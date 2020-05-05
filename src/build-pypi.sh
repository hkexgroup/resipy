#!/bin/bash
# build the python package for pypi

# move out some executable out first
mkdir resipy/exeb
mv resipy/exe/cR2.exe resipy/exeb/cR2.exe
mv resipy/exe/cR3t.exe resipy/exeb/cR3t.exe
mv resipy/exe/R3t.exe resipy/exeb/R3t.exe

# compile source distribution
python3 setup.py sdist bdist_wheel
# only the source distribution allows to run post install script to download binaries

# restore exe
mv resipy/exeb/cR2.exe resipy/exe/cR2.exe
mv resipy/exeb/cR3t.exe resipy/exe/cR3t.exe
mv resipy/exeb/R3t.exe resipy/exe/R3t.exe
rm -r resipy/exeb

#python3 -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*
#python3 -m twine upload dist/*
