#cythonize meshCalc.pyx
from setuptools import setup, Extension 
from Cython.Build import cythonize
import numpy as np
import platform
import os, shutil 

if platform.system() == 'Linux': #open mp support on linux 
    ext_modules = [
        Extension(
            "meshCalc",
            ["meshCalc.pyx"],
        ),
        Extension(
        	"recipCalc",
        	["recipCalc.pyx"]
        	)
    ]
elif platform.system() == 'Windows': #open mp support on windows
    ext_modules = [
        Extension(
            "meshCalc",
            ["meshCalc.pyx"],
        ),
        Extension(
        	"recipCalc",
        	["recipCalc.pyx"]
        	)
    ]
else: # macOS
    '''IMPORTANT
		Assuming "Homebrew" is already installed. To successfully compile meshCalc on macOS
		you need to install gcc through Homebrew by:
		brew install cmake gcc
        brew install cmake gcc@xx for specific version of gcc - then change below to specific version
        Note: gcc-13 works on Apple Silicon
        
        Edit: With the openmp dependency removed, extension can be compiled with apple clang 
    '''
    os.environ['CC'] = 'gcc' # gcc-xx based on gcc version
    os.environ['CXX'] = 'g++' # g++-xx based on g++ version
    ext_modules = [
        Extension(
            "meshCalc",
            ["meshCalcUnix.pyx"],
        ),
        Extension(
        	"recipCalc",
        	["recipCalc.pyx"]
        	)
    ]

setup(
    ext_modules=cythonize(ext_modules),
    include_dirs=[np.get_include()]
)   

# move compiled file into cextension folder 
for f in os.listdir():
    if f.endswith('.pyd') or f.endswith('.so'):
        shutil.move(f, os.path.join('cext',f))
        

#run in console under working directory 
#"python cbuild.py build_ext --inplace"

