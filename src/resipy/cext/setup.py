#cythonize meshCalc.pyx
from distutils.core import setup#, Extension 
from Cython.Build import cythonize
import numpy as np
setup(
    ext_modules=cythonize("meshCalc.pyx"),
    include_dirs=[np.get_include()]
)   

#run in console under working directory 
#"python setup.py build_ext --inplace"

