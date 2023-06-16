# building ResIPy app
# IMPORTANT NOTE: For now pyvista <= 0.28, pyvistaqt <= 0.6 and vtk == 8.1.2 on python 3.7 only work

pyinstaller -y ui-osx.spec
mv ./dist/ResIPy.app ./ResIPy-macos/ResIPy.app
hdiutil create ./dist/ResIPy-macos.dmg -srcfolder ResIPy-macos -ov
mv ./ResIPy-macos/ResIPy.app ./dist/ResIPy-macos.app
