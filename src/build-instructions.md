# Building Instructions
2018-09-30


# 1. Virtual Environment
You need to work inside a python virtual environment. Use the `virtualenv` command to create a new one. Use python 3.X.
Create it using:
``` shell
mkdir <my_virtual_env>
virtualenv <my_virtual_env> -p python3 # use Python 3 if it's not the default on your system
<my_virtual_env>\Scripts\activate.bat #windows version
source <my_virtual_env>\bin\activate #linux version
```
Note: from `virtualenv>12.0.0` there is a bug (https://github.com/pyinstaller/pyinstaller/issues/4064) and `distutils` is not include anymore
in the virtualenv. The work-around is to downgrade virtualenv to 12.0.0 using pip and then create the virtual enviroment.
```
pip install virtualenv==12.0.0
```


# 2. Pyinstaller and packages
See doc online. Activate your virtual environment (as above). Inside this environnement install only the packages needed (see list below) and the *lastest* version of pyinstaller 3.4 using pip. 
`pip install <package_name>`

Packages needed:
- numpy
- matplotlib
- pandas
- pyinstaller
- PyQt5


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




