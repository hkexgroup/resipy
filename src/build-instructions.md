# Building Instructions
Dated: 2019-08-05

# 1. Virtual Environments
You should work inside a python virtual environment in order to build ResIPy into a frozen executable.
There are a few options for initialising a virtual environment, the first is to 
use a the `virtualenv` package, or through the vanilla python `venv` command. Additionally, you could use `conda env` through Anaconda or Miniconda.

## using virtualenv
Use python > 3.6.X. If not installed already then install it through pip*. 

*Note: from `virtualenv>12.0.0` there is a bug (https://github.com/pyinstaller/pyinstaller/issues/4064) and `distutils` is not include anymore
in the virtualenv. The work-around is to downgrade virtualenv to 12.0.0 using pip and then create the virtual enviroment.
```
pip install virtualenv==12.0.0
```
To create and activate a python virtual environment through `virtualenv` : 
``` shell
mkdir <my_virtual_env>
virtualenv <my_virtual_env> -p python3 # use Python 3 if it's not the default on your system
<my_virtual_env>\Scripts\activate.bat #windows version
source <my_virtual_env>\bin\activate #linux version
```
## using python venv command
Using python > 3.6.X. This method can be advantageous as it is results in overall 
smaller compiled executables for windows. 

``` shell
python -m venv <my_virtual_env> # this creates the virtual environment folder in the working directory. 
<my_virtual_env>\Scripts\activate.bat #windows version
source <my_virtual_env>\bin\activate #linux version
```
Note: As we are using pyinstaller to compile resipy we may also need to install additional dependencies on the OS level. 

See: https://pyinstaller.readthedocs.io/en/stable/requirements.html for requirements needed to get going with pyinstaller. 

On linux you'll likley need to have `binutils` installed if not already: 
```shell
sudo pacman -S binutils # if using arch based linux distro 
sudo apt-get install binutils # if using debian/ubuntu based distro 
```

Conda virtual environment: 
```shell
conda create -n <env_name> python
conda activate <env_name>
```

# 2. Pyinstaller and packages
Activate your virtual environment (as above). Inside this environnement install only the packages needed (see list below) and the *lastest* version of pyinstaller (3.5) using pip. 
`pip install <package_name>`

Packages needed:
- numpy
- scipy
- matplotlib
- pandas
- pyinstaller
- PyQt5

Inside the src directory of this reprository there are the following files documenting the package versions required to compile ResIPy:
 - "build_windows_packages.txt" 
 - "build_linux_packages.txt" 

Pay careful attention to the version numbers, as various versions of different packages have known incompatiblities with pyinstaller; particularly the latest version of numpy.
Specific packages and versions can be installed through pip as follows: 
`pip install <package_name>==<version> <package_name>==<version> ...` 


# 3. Building (3 ways)
First change the console directory the src directory of ResIPy: `cd <path_to_resipy>/pyr2/src`

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
Here pyinstaller is initialising using .spec files for different scenarios which we wrote to embed the relevant executable files and supporting directories.  

**notes: 
- We used pyinstaller version 3.5 and,
- PyQt 5.9.x
- For windows there is a bug with pyinstaller and the conda version of pyqt. liklehood is it wont work. We suggest using winpython to build ResIPy on windows.
- if latest version of setuptools causes issue on Windows, downgraded it
- need pyvista==0.23.1




