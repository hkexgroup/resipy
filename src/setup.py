# to build run : python3 setup.py sdist bdist_wheel
import setuptools
# from setuptools.command.develop import develop
# from setuptools.command.install import install
# import os
# import sys
# from resipy.R2 import ResIPy_version

# class PostDevelopCommand(develop):
#     """Pre-installation for development mode."""
#     def run(self):
#         develop.run(self)
#         installBinaries()

# class PostInstallCommand(install):
#     """Pre-installation for installation mode."""
#     def run(self):
#         install.run(self)
#         installBinaries()

# def checkExe(dirname):
#     import requests # doesn't work if 'requests' not already installed
#     exes = ['R2.exe','cR2.exe','R3t.exe','cR3t.exe','gmsh.exe']
#     print(os.listdir(dirname))
#     for exe in exes:
#         fname = os.path.join(dirname, exe)
#         print(fname, 'not found')
#         if os.path.exists(fname) is not True:
#             print('Downloading ' + exe + '...', end='', flush=True)
#             response = requests.get("https://gitlab.com/hkex/pyr2/-/raw/master/src/resipy/exe/" + exe)
#             with open(fname, 'wb') as f:
#                 f.write(response.content)
#             print('done')

# def installBinaries():
#     for p in sys.path:
#         if os.path.isdir(p) and 'resipy' in os.listdir(p):
#             dirname = os.path.join(p, 'build','lib','resipy', 'exe')
#             checkExe(dirname)
#             return
        
with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="resipy",
    version="2.2.0", # CHANGE HERE
    author="HKEx",
    description="API for resistivity and IP inversion/modelling around R2 codes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://gitlab.com/hkex/pyr2",
    packages=setuptools.find_packages(),
    install_requires=['numpy','matplotlib','pandas','scipy','requests'],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    # cmdclass={'install':PostInstallCommand,
    #            'develop':PostDevelopCommand}
)
