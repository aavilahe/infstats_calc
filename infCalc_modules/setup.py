#!/usr/bin/env python
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
	cmdclass = {'build_ext': build_ext},
	ext_modules = [Extension("infCalc_Calcxs", ["infCalc_Calcxs.pyx"],
					extra_compile_args=["-O3"])]
)
