# Building Instructions
2018-09-30


# 1. Virtual Environment
You need to work inside a python virtual environment. Use the `virtualenv` command to create a new one. Use python 3.X.
Create it using:
```command line
cd my_directory
virtualenv <my_virtual_env>
<my_virtual_env>\Scripts\activate.bat #windows version
source <my_virtual_env>\bin\activate #linux version
```



# 2. Pyinstaller and packages
See doc online. Activate your virtual environment (as above). Inside this environnement install only the packages needed (see list below) and the *lastest* version of pyinstaller 3.4 using pip. 
`pip install <package_name>`

Packages needed:
- numpy
- scipy
- matplotlib
- pandas
- pyinstaller


# 3. Buidling (3 ways)
There are 3 ways of building pyR2:
- single zip file
    1. `pyinstaller ui-dir.spec`
- single exe without spash screen
    1. `pyinstaller ui-exe.spec`
- single exe with splash screen
    1. `pyinstaller ui-dir.spec`
    2. manually zip and move the directory `r2gui/dist/ui` to `r2gui/ui.zip`
    3. `pyinstaller splashScreen-exe.spec`

The advantage of the splash screen is that it allows the user to see what is going on when the self-extracting single exe is expanding in a temporary directory. It should prevent too impatient user to click too many times on the app thinking it didn't start.

Also for now there is still the console that is shown when starting the app. But for release it will be hidden.

**notes: 
- we used pyinstaller version 3.4 and,
- pyqt 5.9.x
- for windows there is a bug with pyinstaller and the conda version of pyqt. liklehood is it wont compile. We suggest using winpython to build pyR2 on windows. 




