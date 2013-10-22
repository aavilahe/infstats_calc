infCalc.py
==========

infCalc.py is simple tool that calculates various information-based
statistics for pairs of columns between two protein alignments.


Authors
=========

Aram Avila-Herrera (Aram.Avila-Herrera at ucsf dot edu)


INSTALL
========
For now... just edit the helper script to export `__SRC_PATH`
```bash
#export __SRC_PATH="/path/to/infstats"
export __SRC_PATH="/path/to/where/you/downloaded/infCalc"
```

and call normally
```bash
/path/to/dir/runInfCalc.sh [ args ]

```


Files
=========

runInfCalc.sh -- helper script until I figure out how to package this
infCalc.py -- the main script
infCalc_Modules/infCalc_Aux.py -- extra objects and functions
infCalc_Modules/infCalc_Calcxs.pyx -- extra cython functions

test/ -- directory with example files for a test run


NOTES
========
	parallelized in-mem p-value calculation removed
	under-the-hood functions cythonized


-----------------
Tasks:
-----------------
			-[x] remove sims
			-[x] remove orgdb loading (for sequence grouping)
			-[x] remove parallelized bootstrapping code
			-[x] remove grouping/reweighting
			-[x] add docstrings
			-[x] starting with main.. unbreak what is left
			- make remove_gapped_sites() an option
	        - figure out how to package

