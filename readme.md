GUI for R2 family code
======================

This python wrapper around the R2 family code (http://www.es.lancs.ac.uk/people/amb/Freeware/R2/R2.htm)
provides a standalone graphical user interface (GUI) along with a python API for use in jupyter notebook.


Project structure
-----------------

![project structure](structure.png)*General overlook of the project structure with the three main parts.*

Clone
-----
To clone the project please use git


Examples
--------

using the jupyter notebook:
``` python
    import r2gui.api as r2
    r2.readData(filename)
    mesh = r2.mesh()
    mesh.show() # display mesh using matplotlib
    r2.filtering(minErr=0)
    r2.invert()
    r2.showSection()
```




Best practices
--------------

Here are a set of best coding practices that we need to respect in order to make a
simple, consistent and as easy to maintain as possible code:

- use classes to group functions inside the python API (for instance all functions dealing with meshes should be implemented inside a big Mesh class if possible)
- document functions : you can document your function directly in the code using the ReStructuredText convention (<link needed>) or just use comment with #
- separation of API and GUI should allow to use the python API in jupyter notebook as well as in standalone GUI using pyQt5





