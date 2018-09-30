# Building Instructions
2018-09-30


# 1. Virtual Environment
You need to work inside a python virtual environment. Use the `virtualenv` command to create a new one. Use python 3.X


# 2. Pyinstaller and packages
See doc online. Inside this environnement install only the packages needed (see list below) and the *last development* version of pyinstaller.

Packages needed:
- numpy
- scipy
- matplotlib
- pandas


# 3. Buidling (3 ways)
There are 3 ways of building pyR2:
- single zip file
    1. `pyinstaller ui-dir.spec`
- single exe without spash screen
    1. `pyinstaller ui-exe.spec`
- single exe with splash screen
    1. `pyinstaller ui-dir.spec`
    2. manually zip the directory `dist/ui` to `dist/ui.zip`
    3. `pyinstaller splashScreen-exe.spec`

The advantage of the splash screen is that it allows the user to see what is going on when the self-extracting single exe is expanding in a temporary directory. It should prevent too impatient user to click too many times on the app thinking it didn't start.

Also for now there is still the console that is shown when starting the app. But for release it will be hidden.




