Contributing
===

ResIPy is an open-source software. Therefore we welcome contribution from all developpers.


Best practices
--------------
Here are a set of best coding practices that we need to respect in order to make a
simple, consistent and as easy to maintain as possible code:

- use **classes** to group functions inside the Python API (for instance all functions dealing with meshes should be implemented inside a big Mesh class if possible)
- **document** functions : you can document your function directly in the code using [Numpy docstring](https://numpydoc.readthedocs.io/en/latest/format.html) or just use comment with # for line specific comments
- **separation of API and GUI** should allow to use the Python API in jupyter notebook as well as in standalone GUI using PyQt5


Use of git
----------
Below is the usual commands you are likely to go through if you contribute to this project.

First ensure you have cloned the project and are in the main project directory.
```bash
git clone https://gitlab.com/hkex/resipy
cd resipy
```
Second, you can either (1) create a new branch for your changes (recommended) or (2) use the default `develop` branch.
If you choose (1) you can create a new branch using:
```bash
git checkout -b <name_of_branch>
```

If you choose (2) make sure `develop` is ahead of `master` before changing anything. You can also use those instructions
to update your own branch with `master`.
```bash
git checkout develop # change to develop branch
git fetch origin # get all changes from master
git merge origin/master # merge those change to the develop branch
git push # push them to the online repository
```
 
Finally the typical workflow is as following:
1. Change branch: `git checkout <name_of_branch>`. You can see which branch you use by typing:`git branch`
2. Check you are up to date with the master branch (as shown above): `git fetch origin; git merge origin/master`
3. Operates your changes in the files
4. Use `git status` to see which file need to be added to the repository and which files have just been changed.
5. Use `git add newfile` to add new files.
6. Use `git commit -a` to add a commit messages to the changes you are about to push.
7. Use `git push origin <name_of_branch>` to push your changes to the repository.
8. Go on gitlab and on the project page you will see an invitation to create a merge
9. request from the branch you have just pushed to. You can also go to Repository > 
10. Branches and create a merge request from the `<name_of_branch>` branch.


API structure
-----

This document describes the structre of the API and the files in it. Please refer
to it before doing major modification to the code.


- `src/ui.py`: GUI using PyQt5, it calls `R2 class` for processing
- `src/resipy/Project`: main class that provides import, filtering, mesh creation and inversion,
it calls functions from others scripts for specific purposes. The code is divided in tabs.
- `src/resipy/Survey.py`: important class that contains all data relative to one survey.
It provides methods for looking for reciprocal, fitting error model and filtering data
- `src/resipy/parsers.py`: contains parser function to read data file from different instrument (e.g.: SyscalPro)
- `src/resipy/meshTools.py`: contains the `mesh` class for all mesh creation (quadrilateral and triangular) and plotting
- `src/resipy/gmshWrap.py`: provides a wrapper for `gmsh.exe` to create triangular mesh
- `src/resipy/DCA.py`: provides function for DCA filtering of IP decay curves
- `src/resipy/r2in.py`: provides function to write `R2.in`, `cR2.in`, `R3t` and `cR3t.in` file. Is called by the `Project` class.


Building the software
----

Bundles are build using `pyinstaller` (`pip install pyinstaller`). Different types of bundles can be produced:
- `.zip` bundle, the user unzip it and inside can run `ResIPy.exe` (Windows) or just `ResIPy` (linux), or extractthe `.app` (mac)
- `.exe` (Windows), no extension (linux), `.dmg` (mac). Windows and Linux are self-extractable archives which auto-extract themselves to a temporary directory at run-time (while the splashscreen is loading). Mac dmg image can be mounted and the .app can be drag and dropped to the Applications folder.

Building is automated by the `build.sh` (linux), `build.bat` (windows), `build-app.sh` (mac) scripts.
An additional `build-pypi.sh` is supplied to build PyPi packages.

Note that you need to setup a python environment `pyenv`, as referenced by the scripts, which contains the necessary packages to run ResIPy (see requirements.txt and requirments-gui.txt). `pyenv` should be located at the same level as the repository folder. The scripts should be run from withing the `src` folder.

e.g.
```
pyenv/
resipy/
    doc/
    jupyter-notebook/
    src/
        build.sh
        build.bat
        build-app.sh
        build-pypi.sh
```

To build the Windows bundle, do:
```sh
cd resipy/src
build.bat
```

The bundles produced will be located in `resipy/src/dist/`.



Please see the ![documentation](https://hkex.gitlab.io/resipy/api.html) for a full description of the API.

        
