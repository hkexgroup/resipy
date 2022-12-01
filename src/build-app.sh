# source ../../pyenv/bin/activate


# building ResIPy app
pyinstaller -y ui-osx.spec
mv ./dist/ResIPy.app ./macdmg/ResIPy.app
hdiutil create ./dist/ResIPy-macos.dmg -srcfolder macdmg -ov
mv ./macdmg/ResIPy.app ./dist/ResIPy-macos.app
