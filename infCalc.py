#!/usr/bin/python
''' infCalc.py -- Estimates marginal and joint entropies for
	pairs of alignment columns

	Calculates various information statistics based
	on marginal and joint entropies

	Output is tab delimited

	Cython optimized

	!! Removed in-memory bootstrapping
	!! Removing sequence grouping

'''

import sys
import getopt
from math import isnan

# module path ... fix with egg install later?
__SRC_PATH = '/home/aram/dev-src/infStats_calc/'
sys.path.append(__SRC_PATH+'infCalc_modules')

import infCalc_Aux
import infCalc_Calcxs


# Global Constants
EPS = 2e-6

# input
def getopts(args):
''' Parse command line or config file for options and return a dict.

	using getopt for compatibility

'''

	optlist, args = getopt.getopt(args, 'ho:c:T:', 
							['help', 'outdir=', 'control_file=',
							'vir_aln=', 'host_aln=',
							'vir_sim=', 'host_sim=',
							'vir_keep=', 'host_keep=',
							'taxon_pairs='])

	options = dict()
	usage = 'usage: %s [-h | -c ctl_file | '+\
			'--vir_aln=align1 --host_aln=align2 '+\
			'--vir_sim=sim1 --host_sim=sim2 '+\
			'--vir_keep=keep1 --host_keep=keep2 '+\
			'--taxon_pairs=virushost.pair ] '+\
			'[-o | --outdir directory] '+\
			'[-T num_threads]'

	# read command line
	for opt, val in optlist:
		if opt in ('-h', '--help'):
			sys.exit(usage)
		if opt in ('-o', '--outdir'):
			options['outdir'] = val
		if opt in ('-c', '--control_file'):
			options['ctl'] = val
		if opt == '--vir_aln':
			options['vir_aln'] = val
		if opt == '--host_aln':
			options['host_aln'] = val
		if opt == '--vir_sim':
			options['vir_sim'] = val
		if opt == '--host_sim':
			options['host_sim'] = val
		if opt == '--vir_keep':
			options['vir_keep'] = val
		if opt == '--host_keep':
			options['host_keep'] = val
		if opt == '--taxon_pairs':
			options['taxon_pairs'] = val
		if opt == '-T':
			options['num_threads'] = val
	
	# read control file
	if 'ctl' in options:
		options = parse_ctl(options['ctl'], options)

	# set defaults
	if 'outdir' not in options:
		options['outdir'] = './'
	if 'num_threads' not in options:
		#options['num_threads'] = min(cpu_count() / 2, 1)
		options['num_threads'] = 1 # multiprocessing module NOT loaded
	else:
		options['num_threads'] = int(options['num_threads'])
	
	if not set(['vir_aln', 'host_aln',
			'vir_sim', 'host_sim',
			'vir_keep', 'host_keep', 'taxon_pairs']) <= set(options.keys()):
			sys.exit(usage)
	
	# print options
	print >>sys.stderr, 'running with options:'
	[ sys.stderr.write('\t%s: %s\n'%(opt, val)) 
					for opt, val in sorted(options.iteritems()) ]

	return options

def parse_ctl(ctl_fn, options):
''' Parse config file and return updated options dict


'''

	print >>sys.stderr, 'reading options from: %s'%ctl_fn
	ctl_fh = open(ctl_fn, 'r')
	for line in ctl_fh:
		line = line.strip()
		if line.startswith('#') or '\t' not in line:
			continue
		opt, val = line.split('\t')[:2]
		if opt in options:
			continue
		else:
			options[opt] = val
	
	return options

def read_sites(sites_fn):
''' Reads space delimited file and return list of sites to process.

	Sites should be 0-indexed non-negative integers.

'''

	if sites_fn.lower() == 'all':
		return 'all'
	sites_fh = open(sites_fn, 'r')
	return map(int, ''.join(sites_fh.readlines()).split())

def remove_gapped_sites(keep_sites, aln, taxonPairs):
''' From a list of sites, remove those with gaps and return a list.

'''

	# this patch replaces 'all' keyword with range
	if keep_sites == 'all':
		keep_sites = range(aln.num_cols)

	# remove gapped sites
	remove_these = set()

	# make a set of all taxons
	taxonSet = set()
	[ [ taxonSet.add(x) for x in (v,h)] for (v,h) in taxonPairs ]
	for site_i in keep_sites:
		col = aln.get_site(site_i)
		for orgID in col.keys():
			if '-' in col[orgID] and orgID in taxonSet:
				#print 'removing site', site_i, 'because', orgID
				remove_these.add(site_i)
				break
	return sorted(list(set(keep_sites) - remove_these))

def get_jobname(vir_aln, host_aln):
''' Get jobname from alignment filenames
	
	No-extension basename concatenation

'''

	return vir_aln.rsplit('/', 1)[1].rsplit('.phy', 1)[0] +\
			'_' +\
			host_aln.rsplit('/', 1)[1].rsplit('.phy', 1)[0]

# output
def print_output(resObj, out_fn):
''' Print information stats for each pair of sites in resObj

	resObj is a dict.
		keys are tuples of sites
		values are lists of stats

'''
	outfh = open(out_fn, 'w')
	print >>outfh, 'h1, h2, hj, mi, vi, mi_Nhmin, mi_Nhj'.replace(', ', '\t')
	statstr = '\t%.6f' * 7 
	for coord, stat in resObj.iteritems():
		#print >>outfh, str(coord) + '\t' + '\t'.join(map(str, stat))
		print >>outfh, str(coord) + statstr % ( stat )
	outfh.close()

# meat
def doCalc(vir_aln, host_aln, vir_keep, host_keep, taxonPairs):
	''' loops through columns of alignments
		passes each alignment column labeled with organism name
		to makeProbabilities() as dict k = id, v = symbols

		taxonPairs tells makeProbabilities how to match up
		viruses with hosts

	'''

	if vir_keep == 'all':
		vir_keep = range(vir_aln.num_cols)
	if host_keep == 'all':
		host_keep = range(host_aln.num_cols)
	
	stats = dict()
		
	for i in vir_keep:
		vir_aln_col = vir_aln.get_site(i)
		for j in host_keep:
			host_aln_col = host_aln.get_site(j)
			(m1, m2, jd) = miCalcxs.makeProbabilities(vir_aln_col, host_aln_col, taxonPairs)
			# stats should now have an extra member (VI after MI)
			stats[(i, j)] = miCalcxs.calcStats(m1, m2, jd)
	return stats

def main(options):
	# make output names
	jobname = get_jobname(options['vir_aln'], options['host_aln'])
	pval_fn = options['outdir'] + '/' + jobname + '.pval'
	out_fn = options['outdir'] + '/' + jobname + '.out'

	# load alignments
	vir_aln, vir_orgdb = miAux.read_org_and_phy(options['vir_aln'])
	host_aln, host_orgdb = miAux.read_org_and_phy(options['host_aln'])
	print >>sys.stderr, "alignments loaded"

	# load virus-host pairings map
	taxonPairs = miAux.read_taxonPairs(options['taxon_pairs'], vir_orgdb.values(), host_orgdb.values())
	print >>sys.stderr, "virus-host pairings read"

	# load list of sites to compare
	vir_keep = read_sites(options['vir_keep'])
	host_keep = read_sites(options['host_keep'])
	print >>sys.stderr, "site lists loaded"

	vir_keep = remove_gapped_sites(vir_keep, vir_aln, taxonPairs)
	host_keep = remove_gapped_sites(host_keep, host_aln, taxonPairs)
	print >>sys.stderr, "site lists degapped"
	
	print 'vir_keep (max) = %d; ncol = %d'%( sorted(vir_keep)[-1], vir_aln.num_cols)
	print 'host_keep (max) = %d; ncol = %d'%( sorted(host_keep)[-1], host_aln.num_cols)

	# calculate info stats

	resObj = doCalc(vir_aln, host_aln, vir_keep, host_keep, taxonPairs)
	print >>sys.stderr, "main calculations completed"

	print_output(resObj, out_fn)
	print >>sys.stderr, "done"

if __name__ == "__main__":
	print "DEBUG MODE"
	print sys.argv
	options = getopts(sys.argv[1:])
	main(options)

