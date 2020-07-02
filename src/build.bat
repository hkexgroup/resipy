..\..\pyenv\Scripts\activate.bat && ^
pyinstaller -y ui-dir.spec && ^
mkdir dist\ui\Lib && ^
cd dist && ^
7z a -aoa -r ui.zip ui && ^
cd .. && ^
move /y dist\ui.zip ui.zip && ^
pyinstaller -y splashScreen-exe.spec && ^
move /y dist\ResIPy-launch.exe dist\ResIPy-windows.exe && ^
move /y ui.zip dist\ResIPy-windows.zip
