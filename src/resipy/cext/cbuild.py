from distutils.core import setup, Extension 
import numpy as np
#build 'cythonised c' code with pure python (no need for cython)
#useful for users without cython installed, the module needs to be compiled for each version of python 
#NB: windows users will need visual studio installed to compile (or MS SDK)
#NB: Unix users will need gcc installed 

ext_modules = [
    # Extension(
    #     "meshCalc",
    #     ["meshCalc.c"],
    #     extra_compile_args=['/openmp'],
    #     extra_link_args=['/openmp'],
    # ),
    Extension(
       	"fastRecip",
       	["fastRecip.pyx"]
   	),
]

setup(
      ext_modules=ext_modules,
      include_dirs=[np.get_include()]
)

#run in console under working directory 
#"python cbuild.py build_ext --inplace"