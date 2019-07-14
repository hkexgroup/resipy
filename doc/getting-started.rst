Getting Started
===============

EMagPy aims to process frequency domain electromagnetic measurements (FDEM) taken with electromagnetic induction (EMI) device. EMagPy has a python application programming interface (API) to use in Jupyter-notebook (see section Examples) and a standalone graphical user interface (GUI).

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

