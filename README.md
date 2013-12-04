## infCalc.py ##

infCalc.py is simple tool that calculates various information-based
statistics for pairs of columns between two protein alignments.


## Authors ##

Aram Avila-Herrera (Aram.Avila-Herrera at ucsf dot edu)


## Install ##
Install anywhere you have permission (eg. `~/path/to/infstats_calc`).

1. **Compile the `infCalc_Calcxs` module.**
Requires cython: *<http://docs.cython.org/src/quickstart/install.html>*.
You may skip this step if the included `infCalcxs.so` is already compiled for your system.

	```bash
	cd ~/path/to/infstats_calc/infCalc_modules/
	bash build_modules.sh
	```

2. **Edit the helper script.** Edit `__SRC_PATH` in `runInfCalc.sh` to point to `~/path/to/infCalc`

	```bash
	#export __SRC_PATH="/arams/path/to/infstats"
	export __SRC_PATH="~/path/to/infstats_calc"
	```
	
	Call `runInfCalc.sh` from the install directory or place in your `${PATH}`
	to call normally.	
	
	```bash
	# let's try the test example
	cd ~/path/to/infstats_calc/test_dir/
	../runInfCalc.sh -c test.ctl
	# or after placing runInfCalc.sh in ${PATH}
	runInfCalc.sh -c test.ctl
	```


## Files ##

- runInfCalc.sh -- helper script until I figure out how to package this
- infCalc.py -- the main script
- infCalc_modules/infCalc_Aux.py -- parsing functions (python)
- infCalc_modules/infCalc_Calcxs.pyx -- functions for probabilities and information stats (cython)
- test_dir/ -- directory with example files for a test run

-----------------
##### TODO: #####
- [x] remove sims
- [x] remove orgdb loading (for sequence grouping)
- [x] remove parallelized bootstrapping code
- [x] remove grouping/reweighting
- [x] add docstrings
- [x] starting with main.. unbreak what is left
- [x] make remove_gapped_sites() an option
- [ ] add sequence weighting
- [ ] figure out how to package

