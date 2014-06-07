## infCalc.py ##

infCalc.py is simple script that calculates various information-based
statistics for pairs of columns between two protein alignments.


## Authors ##

Aram Avila-Herrera (Aram.Avila-Herrera at ucsf dot edu)


## Install ##
Install anywhere you have permission (eg. `~/path/to/infCalc`).

1. **Compile the `infCalc_Calcxs` module.**
Requires cython: *<http://docs.cython.org/src/quickstart/install.html>*.
You may skip this step if the included `infCalcxs.so` is already compiled for your system.

	```bash
	cd ~/path/to/infCalc/infCalc_modules/
	bash build_modules.sh
	```

2. **Edit the helper script.** Edit `__SRC_PATH` in `runInfCalc.sh` to point to `~/path/to/infCalc`

	```bash
	#export __SRC_PATH="/arams/path/to/infCalc"
	export __SRC_PATH="~/path/to/infCalc"
	```
	
	Call `runInfCalc.sh` from the install directory or place in your `${PATH}`
	to call normally.	
	
	```bash
	# let's try the test example
	cd ~/path/to/infCalc/test_dir/
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
- [ ] add sequence weighting
- [ ] package properly

