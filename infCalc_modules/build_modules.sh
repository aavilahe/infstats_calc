#!/bin/bash
# runs setup.py to compile cython modules
#
#

# production ready
python setup.py build_ext --inplace

# test optimization
#cython -a `ls ./*.pyx`

