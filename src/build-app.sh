# building ResIPy app

pyinstaller -y ui-osx.spec
mv ./dist/ResIPy.app ./ResIPy-macos/ResIPy.app
hdiutil create ./dist/ResIPy-macos.dmg -srcfolder ResIPy-macos -ov
mv ./ResIPy-macos/ResIPy.app ./dist/ResIPy-macos.app
