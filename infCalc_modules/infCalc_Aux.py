#!/usr/bin/env python
''' infCalc_Aux.py -- Objects and functions for infCalc.py

	class seqAln
	read_phy()
	read_seqID_pairs()

'''

__author__ = 'Aram Avila-Herrera'
__email__ = 'Aram.Avila-Herrera@ucsf.edu'

import sys
import re

class seqAln:
	''' seqAln -- A sequence alignment stored as a list of dicts.

		The alignment is stored as a list of sites represented as dicts whose key is the sequence
		identifier, and the value is the sequence at that site for that sequence.

		eg:
		seq1       ACGGF
		seq2       ACGDF
		seq3       ATGEE

		site[0] := {seq1 : A, seq2: A, seq3: A}
		site[1] := {seq1 : C, seq2: C, seq3: T}
		site[2] := {seq1 : G, seq2: G, seq3: G}
		
	'''

	def __init__(self, num_seqs, num_cols):
		self.num_cols = int(num_cols)
		self.num_seqs = int(num_seqs)
		self.sites = [dict() for i in xrange(self.num_cols)] # list of dicts{key=seqID, symbol}
		#self.orgdb = orgdb # dict{key=seqID; value=orgname}

	def store(self, seqID, seq):
		''' add a sequence to the alignment (one site at a time)

		'''

		for pos, char in enumerate(seq):
			if seqID in self.sites[pos]:
				self.sites[pos][seqID] += char
			else:
				self.sites[pos][seqID] = char

	def get_seqIDs(self):
		return self.sites[0].keys()
	
	def get_site(self, i):
		''' returns a site as a dict{key=seqID, symbol}

		'''

		return self.sites[i]
	
def read_phy(fh, orgdb=None):
	''' reads phylip file and returns seqAln object.

	'''

	seqLine = False
	for line in fh:
		line = line.rstrip('\n\r')
		if line == '':
			continue
		if seqLine:
			seqID = line[:10].strip()
			seq = re.sub('\s+' ,'', line[10:])
			aln.store(seqID, seq)
		else:
			num_seqs, num_cols = line.split()[:2]
			aln = seqAln(num_seqs, num_cols)
			seqLine = True
	return aln

def read_org(fh):
	''' deprecated. -- reads seqID to organismID file

	'''

	orgdb = dict()
	for line in fh:
		line = line.rstrip('\n\r')
		seqid, orgid = line.split('\t')[:2]
		orgdb[seqid.strip()] = orgid.strip() # strip just in case
	return orgdb

def read_org_and_phy(fn):
	''' deprecated. -- reads phylip and orgdb file

	'''

	orgdb_fn = fn.rsplit('.',1)[0] + '.org_db'
	orgdb_fh = open(orgdb_fn, 'r')
	orgdb = read_org(orgdb_fh)
	orgdb_fh.close()
	aln = read_phy(open(fn, 'r'), orgdb)
	return aln, orgdb

def read_sim(fh):
	''' deprecated. -- reads pickled simulations (1000 alignments)

	'''

	import cPickle

	simObj =  cPickle.load(fh)
	if len(simObj) == 0:
		print >>sys.stderr, "sim load failed"

	print >>sys.stderr, "sim read"	
	return simObj

def read_seqID_pairs(fn):
	''' Reads two-column tab-delimited file of paired seqIDs
		and returns a list of tuples of seqIDs.

	'''

	seqID_pairs = list()
	for line in open(fn, 'r'):
		line = line.strip()
		if '\t' in line and not line.startswith('#'):	
			seqID_pairs += [ tuple(line.split('\t')[:2]) ]
	
	return seqID_pairs

def keep_common_seqID_pairs(seqID_pairs, left_seqIDs, right_seqIDs):
	''' Returns list of seqID_pairs in both left_seqIDs and right_seqIDs

	'''
	
	seqIDs_sl = set(left_seqIDs)
	seqIDs_sr = set(right_seqIDs)
	
	common_seqIDs = [ (seqID_l, seqID_r) for (seqID_l, seqID_r) in seqID_pairs
						if seqID_l in seqIDs_sl and seqID_r in seqIDs_sr ]
	return common_seqIDs
	

	

