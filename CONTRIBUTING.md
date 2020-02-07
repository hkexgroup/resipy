Contributing
===

ResIPy is an open-source software. Therefore we welcome contribution from all developpers.


Best practices
--------------
Here are a set of best coding practices that we need to respect in order to make a
simple, consistent and as easy to maintain as possible code:

- use **classes** to group functions inside thePython API (for instance all functions dealing with meshes should be implemented inside a big Mesh class if possible)
- **document** functions : you can document your function directly in the code using the ReStructuredText convention (<link needed>) or just use comment with #
- **separation of API and GUI** should allow to use the Python API in jupyter notebook as well as in standalone GUI using pyQt5


Use of git
----------
Below is the usual commands you are likely to go through if you contribute to this project.

First ensure you have cloned the project and are in the main project directory.
```bash
git clone https://gitlab.com/hkex/pyr2
cd pyr2
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


- `ui.py`: GUI using PyQt5, it calls `R2 class` for processing
- `api/R2`: main class that provides import, filtering, mesh creation and inversion,
it calls functions from others scripts for specific purposes
- `api/Survey.py`: important class that contains all data relative to one survey.
It provides methods for looking for reciprocal, fitting error model and filtering data
- `api/parsers.py`: contains parser function to read data file from different instrument (e.g.: SyscalPro)
- `api/meshTools.py`: contains the `mesh` class for all mesh creation (quadrilateral and triangular) and plotting
- `api/gmshWrap.py`: provides a wrapper for `gmsh.exe` to create triangular mesh
- `api/DCA.py`: provides function for DCA filtering of IP decay curves
- `api/r2in.py`: provides function to write `R2.in` and `cR2.in` file. Is called by the `R2 class`.


Please see the ![documentation](https://hkex.gitlab.io/pyr2/api.html) for a full description of the API.

        
