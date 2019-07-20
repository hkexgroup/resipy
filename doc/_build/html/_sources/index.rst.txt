.. ResIPy documentation master file, created by
   sphinx-quickstart on Wed Aug 29 00:30:45 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

    
ResIPy python API and standalone GUI
====================================

Installation
------------

Clone the gitlab repository::

    git clone https://gitlab.com/hkex/pyr2

To start the GUI from source, navigate through the `src` directory and run `ui.py`::

    cd pyr2/src
    python ui.py
    
From the same `src` directory you can import the module from source using python. Another solution is to install the module from pypi using pip::

    pip install resipy
    
.. note::
    Mac and Linux user will need *wine* to run the inversions.


.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   getting-started
   gui
   api
   auto_examples/index
   


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

   

