source ../../pyenv/bin/activate

# building the cython extension
cd resipy/cext
python setup.py build_ext --inplace
cd ../..

# building ResIPy app
pyinstaller -y ui-osx.spec
mv ./dist/ResIPy.app ./macdmg/ResIPy.app
hdiutil create ./dist/ResIPy-macos.dmg -srcfolder macdmg -ov
mv ./macdmg/ResIPy.app ./dist/ResIPy-macos.app
