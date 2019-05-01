#!/bin/bash
# build the python package as well
python3 setup.py sdist bdist_wheel
#python3 -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*
#python3 -m twine upload dist/*

