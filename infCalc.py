#!/usr/bin/env python
''' infCalc.py -- Calculates various information statistics

	infCalc.py calculates various information statistics based
	on marginal and joint entropies of pairs of protein
	alignment columns.

	Input: Two phylip alignmnents. Control file (optional).
	Output: Tab delimited file of interprotein alignment columns
			and their scores.

'''

__author__ = 'Aram Avila-Herrera'
__email__ = 'Aram.Avila-Herrera@ucsf.edu'



import sys
from os import getenv
from os.path import basename
import getopt
from math import isnan

# module path ... not sure how to fix this yet
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
				'          [-o | --outdir directory]\n'#+\
#				'          [-T num_threads]'
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
	# Make output names
	jobname = get_jobname(options['vir_aln'], options['host_aln'])
	out_fn = options['outdir'] + '/' + jobname + '.out'

	(vir_aln, host_aln,
	vir_keep, host_keep,
	seqID_pairs) = iC_A.load_all_input(options)

	#DBG info
	print "DBG: last column in vir_keep = [ %d ]" % max(vir_keep)
	print "DBG: last column in host_keep = [ %d ]" % max(host_keep)
	print "DBG: num columns in vir_aln = [ %d ]" % vir_aln.num_cols
	print "DBG: num columns in host_aln = [ %d ]" % host_aln.num_cols
	print "DBG: num columns in vir_keep = [ %d ]" % len(vir_keep)
	print "DBG: num columns in host_keep = [ %d ]" % len(host_keep)

	# Calculate info stats
	print "Main calculations...",
	resObj = doCalc(vir_aln, host_aln, vir_keep, host_keep, seqID_pairs)
	print "Completed!"

	print_output(resObj, out_fn)
	print "Done!"

if __name__ == "__main__":
	print >>sys.stderr, "DBG: sys.argv:", sys.argv
	options = getopts(sys.argv[1:])
	main(options)

