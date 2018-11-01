..\pyenv\Scripts\activate.bat
pyinstaller ui-dir.spec
cd dist
7z a -r ui.zip ui
cd ..
move dist\ui.zip ui.zip
pyinstaller splashScreen-exe.spec