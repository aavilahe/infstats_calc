#!/usr/bin/python
''' alignment object
	
	and other methods...???

'''

import cPickle
import sys

class seqAln:
	def __init__(self, num_seqs, num_cols):#, orgdb):
		self.num_cols = int(num_cols)
		self.num_seqs = int(num_seqs)
		self.sites = [dict() for i in xrange(self.num_cols)] # list of dicts{key=seq_ident, value=subcol}
		#self.orgdb = orgdb # dict, key = seqid, value=orgname

	def store(self, seq_ident, seq):
		''' add sequence to alignment

		'''

		for pos, char in enumerate(seq):
			if seq_ident in self.sites[pos]:
				self.sites[pos][seq_ident] += char
			else:
				self.sites[pos][seq_ident] = char

	def get_seqIDs(self):
		return self.sites[0].keys()
	
	def get_site(self, i):
		return self.sites[i]
	
def read_phy(fh, orgdb=None):
	seqLine = False
	for line in fh:
		line = line.rstrip('\n\r')
		if line == '':
			continue
		if seqLine:
			seq_ident = line[:10].strip() # strip just in case
			seq = line[10:].replace(' ', '')
			aln.store(seq_ident, seq)
		else:
			num_seqs, num_cols = line.split()[:2]
			aln = seqAln(num_seqs, num_cols)
			seqLine = True
	return aln

def read_org(fh):
	orgdb = dict()
	for line in fh:
		line = line.rstrip('\n\r')
		seqid, orgid = line.split('\t')[:2]
		orgdb[seqid.strip()] = orgid.strip() # strip just in case
	return orgdb

def read_org_and_phy(fn):
	orgdb_fn = fn.rsplit('.',1)[0] + '.org_db'
	orgdb_fh = open(orgdb_fn, 'r')
	orgdb = read_org(orgdb_fh)
	orgdb_fh.close()
	aln = read_phy(open(fn, 'r'), orgdb)
	return aln, orgdb

def read_sim(fh):
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

