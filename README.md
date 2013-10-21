Clean version of information stats calculator
=============================================

Features:
--------
	Modular: parallelized in-mem p-value calculation removed
	Fast: Cython optimization
	Output: tab delimited
	Filters gapped sites:
	**NO** reweighting / grouping of sequences

-----------------
Top down refactor
-----------------
	*DONE*: starting with main..
			. remove sims
			. remove orgdb loading (for sequence grouping)
			. remove parallelized bootstrapping code
			. remove grouping/reweighting
			. add docstrings
			. unbreak what is left
	*TODO*: . make remove_gapped_sites() an option
	        . figure out how to package

