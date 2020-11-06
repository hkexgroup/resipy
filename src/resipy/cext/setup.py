#cythonize meshCalc.pyx
from distutils.core import setup, Extension 
from Cython.Build import cythonize
import numpy as np
import platform
import os

if platform.system() == 'Linux': #open mp support on linux 
    ext_modules = [
        Extension(
            "meshCalc",
            ["meshCalc.pyx"],
            extra_compile_args=["-fopenmp"],
            extra_link_args=['-fopenmp'],
        )
    ]
elif platform.system() == 'Windows': #open mp support on windows
    ext_modules = [
        Extension(
            "meshCalc",
            ["meshCalc.pyx"],
            extra_compile_args=['/openmp'],
            #extra_link_args=['/openmp'],
        )
    ]
else: # macOS
    '''IMPORTANT
		Assuming "Homebrew" is already installed. To successfully compile meshCalc on macOS
		you need to install gcc through Homebrew by:
		brew install cmake gcc
    '''
    os.environ['CC'] = 'gcc-10' # gcc-xx based on gcc version
    os.environ['CXX'] = 'g++-10' # g++-xx based on g++ version
    ext_modules = [
        Extension(
            "meshCalc",
            ["meshCalc.pyx"],
            extra_compile_args=["-fopenmp"],
            extra_link_args=['-fopenmp'],
        )
    ]

setup(
    ext_modules=cythonize(ext_modules),
    include_dirs=[np.get_include()]
)   

#run in console under working directory 
#"python setup.py build_ext --inplace"

