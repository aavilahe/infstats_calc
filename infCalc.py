#!/usr/bin/env python
''' infCalc.py -- Calculates various information statistics

	infCalc.py calculates various information statistics based
	on marginal and joint entropies of pairs of protein
	alignment columns.

	Input: Two phylip alignmnents. Control file (optional).
	Output: Tab delimited file of interprotein alignment columns
			and their scores.

'''

import sys
from os import getenv
from os.path import basename
import getopt
from math import isnan

# module path ... fix with egg install later?
#__SRC_PATH = '/home/aram/dev-src/infStats_calc/'
__SRC_PATH = getenv("__SRC_PATH")
sys.path.append(__SRC_PATH+'/'+'infCalc_modules')
print "DBG:", "__SRC_PATH:", __SRC_PATH

import infCalc_Aux as iC_A
import infCalc_Calcxs as iC_C


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
#							'vir_sim=', 'host_sim=',
							'vir_keep=', 'host_keep=',
							'seqID_pairs='])

	options = dict()
	usage = (
				'usage: %s [-h | -c ctl_file |\n'+\
				'            --vir_aln=align1 --host_aln=align2\n'+\
#				'            --vir_sim=sim1 --host_sim=sim2\n'+\
				'            --vir_keep=keep1 --host_keep=keep2\n'+\
				'            --seqID_pairs=virushost.pair ]\n'+\
				'          [-o | --outdir directory]\n'+\
				'          [-T num_threads]'
			)% sys.argv[0]

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
#		if opt == '--vir_sim':
#			options['vir_sim'] = val
#		if opt == '--host_sim':
#			options['host_sim'] = val
		if opt == '--vir_keep':
			options['vir_keep'] = val
		if opt == '--host_keep':
			options['host_keep'] = val
		if opt == '--seqID_pairs':
			options['seqID_pairs'] = val
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
	if 'host_keep' not in options:
		options['host_keep'] = 'all'
	if 'vir_keep' not in options:
		options['vir_keep'] = 'all'
	
	if not set(['vir_aln', 'host_aln',
#			'vir_sim', 'host_sim',
			'vir_keep', 'host_keep', 'seqID_pairs']) <= set(options.keys()):
			sys.exit(usage)
	
	# print options
	print 'running with options:'
	[ sys.stdout.write('\t%s: %s\n'%(opt, val)) 
					for opt, val in sorted(options.iteritems()) ]

	return options

def parse_ctl(ctl_fn, options):
	''' Parse config file and return updated options dict
	
	
	'''

	print 'reading options from: %s'%ctl_fn
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

def remove_gapped_sites(keep_sites, aln, seqIDs, gap_threshold=0.3):
	''' From a list of sites, remove those with gaps above a
		threshold and return a list.

	'''

	# this patch replaces 'all' keyword with site range
	# this will be moved to read_sites()
	if keep_sites == 'all':
		keep_sites = range(aln.num_cols)

	# remove gapped sites
	remove_these = set()

	# for each aln, remove sites with more than X% gap in relevant seqs
	for site_i in keep_sites:
		col = aln.get_site(site_i)
		filCol = dict([ (seqID, col[seqID]) for seqID in seqIDs ])
		ngaps = ''.join(filCol.values()).count('-')
		nseqs = float(len(filCol))
		if ngaps / nseqs >= gap_threshold:
			remove_these.add(site_i)
	return sorted(list(set(keep_sites) - remove_these))

def get_jobname(vir_aln, host_aln):
	''' Get jobname from alignment filenames
	
	No-extension basename concatenation

	'''

	return basename(vir_aln).rsplit('.phy', 1)[0] +\
			'_' +\
			basename(host_aln).rsplit('.phy', 1)[0]

# output
def print_output(resObj, out_fn):
	''' Print information stats for each pair of sites in resObj

	resObj is a dict.
		keys are tuples of sites
		values are lists of stats

	'''
	outfh = open(out_fn, 'w')
	print >>outfh, '\t'.join(('Virus_Column', 'Mammal_Column',
								'Vir_Entropy', 'Mam_Entropy', 'Joint_Entropy',
								'MutInf', 'VarInf', 'Zmin_MutInf', 'Zjoint_MutInf'))
	statstr = '\t%.6f' * 7
	for coord, stat in resObj.iteritems():
		print >>outfh, '\t'.join(map(str, coord)) + statstr % ( stat )
	outfh.close()

# meat
def doCalc(vir_aln, host_aln, vir_keep, host_keep, seqID_pairs):
	''' Loops through pairs of alignment columns and returns stats in a dict.
			Passes each alignment column labeled with organism name
				to makeProbabilities()
			seqID_pairs tells makeProbabilities() how to match sequences

		stats is a dict.
			keys are tuples of sites
			values are list of stats: [ MutInf, VarInf, ZminH, ZHj ]

	'''

	stats = dict()
		
	for i in vir_keep:
		vir_aln_col = vir_aln.get_site(i)
		for j in host_keep:
			host_aln_col = host_aln.get_site(j)
			(vir_prob, host_prob, joint_prob) = iC_C.makeProbabilities(vir_aln_col, host_aln_col, seqID_pairs)
			stats[(i, j)] = iC_C.calcStats(vir_prob, host_prob, joint_prob)
	return stats

def main(options):
	# make output names
	jobname = get_jobname(options['vir_aln'], options['host_aln'])
	out_fn = options['outdir'] + '/' + jobname + '.out'

	# load alignments
	#vir_aln, vir_orgdb = miAux.read_org_and_phy(options['vir_aln'])
	#host_aln, host_orgdb = miAux.read_org_and_phy(options['host_aln'])
	vir_aln = iC_A.read_phy(open(options['vir_aln'],'r'))
	host_aln = iC_A.read_phy(open(options['host_aln'],'r'))
	print "alignments loaded"

	## load virus-host pairings map
	# load sequence pairings
	seqID_pairs = iC_A.read_seqID_pairs(options['seqID_pairs'])
	print "virus-host pairings read"

	# load list of sites to compare
	vir_keep = read_sites(options['vir_keep'])
	host_keep = read_sites(options['host_keep'])
	print "site lists loaded"

	print >>sys.stderr, "DBG: remove_gapped_sites() necessary for 'all' keyword"
	vir_keep = remove_gapped_sites(vir_keep, vir_aln,
						zip(*seqID_pairs)[0])
	host_keep = remove_gapped_sites(host_keep, host_aln,
						zip(*seqID_pairs)[1])
	print "site lists degapped"

	#DBG info
	print 'vir_keep (max) = %d; ncol = %d'%( sorted(vir_keep)[-1], vir_aln.num_cols)
	print 'host_keep (max) = %d; ncol = %d'%( sorted(host_keep)[-1], host_aln.num_cols)

	# calculate info stats

	resObj = doCalc(vir_aln, host_aln, vir_keep, host_keep, seqID_pairs)
	print "main calculations completed"

	print_output(resObj, out_fn)
	print "done"

if __name__ == "__main__":
	print >>sys.stderr, "DBG: sys.argv:", sys.argv
	options = getopts(sys.argv[1:])
	main(options)

