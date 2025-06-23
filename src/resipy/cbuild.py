#cythonize meshCalc.pyx
from setuptools import setup, Extension 
from Cython.Build import cythonize
import numpy as np
import platform
import os, shutil, platform

# sort out long integer data types for each platform 
if platform.system() == "Windows":
    print("Platform is windows -> changing long data type")
    long_data_type = "long long"
else: 
    long_data_type = "long"

files2compile = ['meshCalc.pyx', 'recipCalc.pyx']
for f in files2compile: 
    fh = open(f, 'r')
    filecontents = fh.read()
    fh.close() 
    tmpfilecontents = filecontents.format(long_data_type)
    fh = open("_"+f, "w")
    fh.write(tmpfilecontents)
    fh.close() 

if platform.system() == 'Linux': #open mp support on linux 
    ext_modules = [
        Extension(
            "meshCalc",
            ["_meshCalc.pyx"],
        ),
        Extension(
        	"recipCalc",
        	["_recipCalc.pyx"],
        	)
    ]
elif platform.system() == 'Windows': #open mp support on windows
    ext_modules = [
        Extension(
            "meshCalc",
            ["_meshCalc.pyx"],
        ),
        Extension(
        	"recipCalc",
        	["_recipCalc.pyx"],
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
            ["_meshCalc.pyx"],
        ),
        Extension(
        	"recipCalc",
        	["_recipCalc.pyx"]
        	)
    ]

setup(
    ext_modules=cythonize(ext_modules),
    include_dirs=[np.get_include()]
)   

# move compiled file into cextension folder and clean up  
for f in os.listdir():
    if f.endswith('.pyd') or f.endswith('.so'):
        shutil.move(f, os.path.join('cext',f))

for f in files2compile: 
    os.remove("_" + f.replace(".pyx",".c"))
    os.remove("_" + f)

#run in console under working directory 
#"python cbuild.py build_ext --inplace"

