..\..\pyenv\Scripts\activate.bat && ^
pyinstaller -y ui-dir.spec && ^
cd dist && ^
7z ao -r ui.zip ui && ^
cd .. && ^
move dist\ui.zip ui.zip && ^
pyinstaller -y splashScreen-exe.spec
