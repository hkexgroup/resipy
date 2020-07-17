#cythonize meshCalc.pyx
from distutils.core import setup, Extension 
from Cython.Build import cythonize
import numpy as np
import platorm

if platorm.system() == 'Linux':
    earg = '-fopenmp'
elif platorm.system() == 'Windows':
    earg = '/openmp'
else:
    earg = ''
    
ext_modules = [
    Extension(
        "meshCalc",
        ["meshCalc.pyx"],
        extra_compile_args=['-fopenmp'],
        extra_link_args=['-fopenmp'],
    )
]
setup(
    ext_modules=cythonize(ext_modules),
    include_dirs=[np.get_include()]
)   

#run in console under working directory 
#"python setup.py build_ext --inplace"

