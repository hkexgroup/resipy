from distutils.core import setup, Extension 
import numpy as np
#build 'cythonised c' code with pure python (no need for cython)
#useful for users without cython installed, the module needs to be compiled for each version of python 
#NB: windows users will need visual studio installed to compile (or MS SDK)
#NB: Unix users will need gcc installed 

setup(
      ext_modules=[Extension('meshCalc',['meshCalc.c'])],
      include_dirs=[np.get_include()]
)

#run in console under working directory 
#"python cbuild.py build_ext --inplace"