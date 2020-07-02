#!/bin/bash
. ../../pyenv/bin/activate # activate virtual environment

# building the cython extension
cd resipy/cext
python setup.py build_ext --inplace
cd ../..

# compiling one folder and zipping it
pyinstaller -y ui-dir.spec
cd dist
zip -r ui.zip ui
cd ..
mv dist/ui.zip ui.zip

# then compile splashcreen that will unzip the zip
pyinstaller -y splashScreen-exe.spec

mv dist/ResIPy-launch dist/ResIPy-linux
mv ui.zip dist/ResIPy-linux.zip

