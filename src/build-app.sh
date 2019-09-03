source ../../pyenv/bin/activate
pyinstaller -y ui-osx.spec
mv ./dist/ResIPy.app ./macdmg/ResIPy.app
hdiutil create ./dist/ResIPy.dmg -srcfolder macdmg -ov
mv ./macdmg/ResIPy.app ./dist/ResIPy.app
