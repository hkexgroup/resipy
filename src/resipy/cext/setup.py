#cythonize meshCalc.pyx
from distutils.core import setup, Extension 
from Cython.Build import cythonize
import numpy as np
import platform

if platform.system() == 'Linux': #open mp support on linux only 
    ext_modules = [
        Extension(
            "meshCalc",
            ["meshCalc.pyx"],
            extra_compile_args=['-fopenmp'],
            extra_link_args=['-fopenmp'],
        )
    ]
else:
    ext_modules = [
        Extension(
            "meshCalc",
            ["meshCalc.pyx"],
        )
    ]

setup(
    ext_modules=cythonize(ext_modules),
    include_dirs=[np.get_include()]
)   

#run in console under working directory 
#"python setup.py build_ext --inplace"

