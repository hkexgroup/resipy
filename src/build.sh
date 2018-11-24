#!/bin/bash
. ../../pyenv/bin/activate
pyinstaller -y ui-dir.spec
cd dist
zip -r ui.zip ui
cd ..
mv dist/ui.zip ui.zip
pyinstaller -y splashScreen-exe.spec
