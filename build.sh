#!/bin/bash
source ../pyenv/bin/activate
pyinstaller ui-dir.spec
cd dist
zip -r ui.zip ui
cd ..
mv dist/ui.zip ui.zip
pyinstaller splashScreen-exe.spec

