#!/bin/bash
activate()
{
    . ../pyenv/bin/activate
}

activate
pip list
pyinstaller --help
pyinstaller ui-dir.spec
