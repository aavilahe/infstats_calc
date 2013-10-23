infCalc.py
==========

infCalc.py is simple tool that calculates various information-based
statistics for pairs of columns between two protein alignments.


Authors
=========

Aram Avila-Herrera (Aram.Avila-Herrera at ucsf dot edu)


INSTALL
========
If the cython module is built correctly for your architecture 
(`x86_64`), skip **step 1**.

1. Build cythonized module 
Requires Cython <http://docs.cython.org/src/quickstart/install.html>
```bash
cd /my/path/to/infCalc/infCalc_modules/
bash build_modules.sh
```

2. For now... just edit the helper script `runInfCalc.sh` to export `__SRC_PATH`
```bash
#export __SRC_PATH="/path/to/infstats"
export __SRC_PATH="/my/path/to/infCalc"
```
and call normally, for example:
```bash
cd /my/path/to/infCalc/test_dir/
../runInfCalc.sh -c test.ctl
```


Files
=========

- runInfCalc.sh -- helper script until I figure out how to package this
- infCalc.py -- the main script
- infCalc_modules/infCalc_Aux.py -- extra objects and functions
- infCalc_modules/infCalc_Calcxs.pyx -- extra cython functions
- test_dir/ -- directory with example files for a test run


NOTES
========
- parallelized in-mem p-value calculation removed
- under-the-hood functions cythonized


-----------------
Tasks:
-----------------
- [x] remove sims
- [x] remove orgdb loading (for sequence grouping)
- [x] remove parallelized bootstrapping code
- [x] remove grouping/reweighting
- [x] add docstrings
- [x] starting with main.. unbreak what is left
- [ ] make remove_gapped_sites() an option
- [ ] figure out how to package

