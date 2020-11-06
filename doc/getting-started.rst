Getting Started
===============

ResIPy aims to process geoeletrical data (DC and IP). It provides a python application programming interface (API) and a standalone graphical user interface (GUI). Resipy provides high-level filtering, error modelling, inversion/forward modelling and post-processing tools.

.. _guiGif:
.. figure:: ../src/image/teaser.gif
    :alt: Animated gif of the GUI in action
    
    Animation of the workflow of creating synthetic data from forward modelling and then inverting them in the GUI.


The same processing can be achieved with the python API::

    from resipy import Project
    k = Project(typ='R2')
    k.createSurvey('examples/dc-2d/syscal.csv')
    k.invert() # invert measurements
    k.showResults() # display inverted pseudo-section

More examples are available in the Gallery of examples.


Installation
------------

The easiest way is to download one of our standalone executable from gitlab (https://gitlab.com/hkex/resipy).


Clone the gitlab repository and run from source::

    git clone https://gitlab.com/hkex/resipy
    cd resipy/src
    python ui.py
    
Alternatively you can install the API part of the module (so no GUI) from pypi using pip::

    pip install resipy
    
    
.. note::
    Mac and Linux user will need *wine* to run the inversions.


Citing ResIPy
-------------
If you use ResIPy for you work, please cite this paper as:

Blanchy G., Saneiyan S., Boyd J., McLachlan P. and Binley A. 2020.
“ResIPy, an Intuitive Open Source Software for Complex  Geoelectrical Inversion/Modeling.”
Computers & Geosciences, February, 104423. https://doi.org/10.1016/j.cageo.2020.104423.


BibTex code::

    @article{blanchy_resipy_2020,
	    title = {{ResIPy}, an intuitive open source software for complex geoelectrical inversion/modeling},
	    issn = {0098-3004},
	    url = {http://www.sciencedirect.com/science/article/pii/S0098300419308192},
	    doi = {10.1016/j.cageo.2020.104423},
	    pages = {104423},
	    journaltitle = {Computers \& Geosciences},
	    author = {Blanchy, Guillaume and Saneiyan, Sina and Boyd, James and {McLachlan}, Paul and Binley, Andrew},
	    urldate = {2020-02-07},
	    date = {2020-02-04},
	    langid = {english}
    }


