Clean version of information stats calculator
=============================================

Features
--------
	Modular: no auto p-values
	Fast: Cython optimization
	Output: tab delimited
	Filters gapped sites:
	NO reweighting / grouping of sequences

-----------------
Top down refactor
-----------------
	*DONE*: starting with main..
			. remove sims
			. remove orgdb
			. remove bootstrapping
			. remove grouping/reweighting
			. add docstrings
			. unbreak what is left
	*TODO*: . make remove_gapped_sites() an option
	        . figure out how to package

