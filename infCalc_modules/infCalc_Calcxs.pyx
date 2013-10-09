#!/usr/bin/python
''' miCalcxs.pyx -- functions for estimating probability distributions,
					entropy, and mutual information

'''

import sys
from itertools import product

cdef extern from "math.h":
	double log(double)
	double fmin(double, double)
	double fmax(double, double)

cdef double eps = 2e-6
cdef double cdef_float_nan = float('nan')

cpdef tuple makeProbabilities(dict column_1, dict column_2, list seqID_pairs):
	''' Calculates marginal probabilities for symbols
		in each column based on observed frequencies
		
		Calculates joint distribution for pairs of symbols
		in both columns

		Force 1-to-1 relationship between sequences in different columns

		seqID_pairs is a list of tuples containing seqIDs for first and second column
		for each seqID_pair, the subcolumn is extracted and the prob distros are updated

		column_1 and column_2 are dicts.
			keys are seqIDs
			values are string of symbols (eg. amino acids ACDEFG...Y)

	'''

	marginal_1 = dict()
	marginal_2 = dict()
	joint = dict()
	
	cdef int num_seqIDs = len(seqID_pairs)

	cdef int seqID_idx
	for seqID_idx in range(num_seqIDs):
		seqID_1, seqID_2 = seqID_pairs[seqID_idx]
		symbols_1 = column_1[seqID_on1] # str
		symbols_2 = column_2[seqID_on2] # str
		#num_obs_1 = len(symbols_1)
		#num_obs_2 = len(symbols_2)

		#updateMarginal(num_obs_1 * num_seqIDs, symbols_1, marginal_1)
		updateMarginal(num_seqIDs, symbols_1, marginal_1)
		#updateMarginal(num_obs_2 * num_seqIDs, symbols_2, marginal_2)
		updateMarginal(num_seqIDs, symbols_2, marginal_2)
		#updateJoint(num_obs_1 * num_obs_2 * num_seqIDs, symbols_1, symbols_2, joint)
		updateJoint(num_seqIDs, symbols_1, symbols_2, joint)

	sanityCheckDist(marginal_1)
	sanityCheckDist(marginal_2)
	sanityCheckDist(joint)

	return (marginal_1, marginal_2, joint)

cdef int updateMarginal(double possibleObs, bytes seqID_Syms, dict marginal):
	''' increment each symbol's probabilities

	'''

	for sym in seqID_Syms:
		exp_sym, weight = expand(sym)
		inc = weight / possibleObs
		for termsym in exp_sym:
			if termsym in marginal:
				marginal[termsym] += inc
			else:
				marginal[termsym] = inc

cdef int updateJoint(double possibleObs, bytes seqID_Syms_1, bytes seqID_Syms_2, dict joint):
	''' increment each pair of symbol's probabilities

	'''

	for sympair in product(seqID_Syms_1, seqID_Syms_2):
		exp_sympair, weight = expand2(sympair)
		inc = weight / possibleObs
		for termpair in exp_sympair:
			if termpair in joint:
				joint[termpair] += inc
			else:
				joint[termpair] = inc

cdef inline tuple expand(bytes nonterm):
	if nonterm == 'X':
		return 'ACDEFGHIKLMNPQRSTVWY', 0.05
	if nonterm == 'B':
		return 'DN', 0.5
	if nonterm == 'Z':
		return 'EQ', 0.5
	return nonterm, 1.0

cdef inline tuple expand2(tuple nonterms):
	e1, w1 = expand(nonterms[0])
	e2, w2 = expand(nonterms[1])
	return product(e1,e2), w1 * w2
	
cpdef int sanityCheckDist(dict probDist):
	''' make sure these probabilities are ok

	'''

	cdef double sumProbs

	for sym, prob in probDist.iteritems():
		probDist[sym] = fmax(0.0, fmin(1.0, prob))
	sumProbs = sum(probDist.values())
	if sumProbs - eps> 1.0:
		sys.stderr.write(str(probDist) + "\n")
		sys.stderr.write(str(sumProbs) + '\n')
		sys.exit("sanityCheckDist failed:: probs > 1")
	if sumProbs < 0.0:
		sys.stderr.write(str(probDist) + "\n")
		sys.exit("sanityCheckDist failed:: probs < 0")

cpdef tuple calcStats(dict marginal_1, dict marginal_2, dict joint):
	''' calcStats calculates
		- entropy 1 (h1)
		- entropy 2 (h2)
		- entropy joint (hj)
		- mutual information (mi)
		- mi / min(h1, h2)
		- mi / hj

	'''
	
	cdef double h1 = calcH(marginal_1)
	cdef double h2 = calcH(marginal_2)
	cdef double hj = calcH(joint)
	cdef double mi = fmax(h1 + h2 - hj, 0.0)
	cdef double vi = fmax(hj - mi, 0.0)
	cdef double mi_Nhmin
	cdef double mi_Nhj
	
	# normalized mis
	cdef double hmin = fmin(h1, h2)
	if hmin > eps:
		mi_Nhmin = mi / hmin
	else:
		mi_Nhmin = cdef_float_nan
	
	if hj > eps:
		mi_Nhj = mi / hj
	else:
		mi_Nhj = cdef_float_nan
	
	return (h1, h2, hj, mi, vi, mi_Nhmin, mi_Nhj)

cdef inline double calcH(dict probDist):
	''' calculate entropy of probDist

	'''

	# probDist is a dict of syms

	cdef double h = 0.0
	cdef double p
	for p in probDist.itervalues():
		if p > eps:
			h -= p * log(p)
	h = fmax(0.0, h)

	return h

cdef int printP(bytes outname, dict pvals):
	outf = open(outname, 'w')
	for c, plist in pvals.iteritems():
		outf.write(str(c) + '\t' + '\t'.join(map(str,plist)) + '\n')
	return 0

cpdef int Pcount(object queue, int numReps, bytes outname):
	cdef int completed = 0
	cdef double pval_inc = 1.0 / numReps
	pvals = dict()
	while True:
		if completed >= numReps:
			return printP(outname, pvals)
		(x, y, pi) = queue.get()
		if (x, y) == ('X', 'X'):
			completed += 1
			continue
		if (x, y) in pvals:
			pvals[(x,y)][pi] += pval_inc
		else:
			pvals[(x,y)] = [0.0, 0.0, 0.0, 0.0] # mi, vi, mi_Nhmin, mi_Nhj
			pvals[(x,y)][pi] += pval_inc
	
#cdef inline double fmax(double a, double b):
#	if a > b:
#		return a
#	return b
#
#cdef inline double fmin(double a, double b):
#	if a < b:
#		return a
#	return b

