source ../../pyenv/bin/activate
pyinstaller ui-osx.spec
mv ./dist/pyR2.app ./macdmg/pyR2.app
hdiutil create ./dist/pyR2.dmg -srcfolder macdmg -ov
mv ./macdmg/pyR2.app ./dist/pyR2.app
