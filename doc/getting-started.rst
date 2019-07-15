Getting Started
===============

ResIPy aims to process geoeletrical data (DC and IP). It provides a python application programming interface (API) and a standalone graphical user interface (GUI). Resipy provides high-level filtering, error modelling, inversion/forward modelling and post-processing tools.

.. _guiGif:
.. figure:: ../src/image/teaser.gif
    :alt: Animated gif of the GUI in action
    
    Animation of the workflow of creating synthetic data from forward modelling and then inverting them in the GUI.


The same processing can be achieved with the python API::

    from resipy import R2
    k = R2()
    k.createSurvey('test/syscalFile.csv')
    k.invert() # invert measurements
    k.showResults() # display inverted pseudo-section

More examples are available in the Example pseudo-section.

