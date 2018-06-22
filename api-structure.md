API structure
=============

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


        
        
