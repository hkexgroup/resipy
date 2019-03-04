#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 17 14:16:42 2019

@author: jimmy
"""
import numpy as np
import sys
sys.path.append(r'/home/jimmy/phd/pyr2/src')

from api.R2 import R2


#%% 
print('-------------Testing Forward DC Modelling ------------')
k = R2(typ='R2')

k.elec = np.c_[np.linspace(0,24, 24), np.zeros((24, 2))]
k.createMesh(typ='quad')
#
#%% full API function
k.addRegion(np.array([[5,-1.5],[10,-1.5],[10,-3.5],[5,-3.5]]), 10, -3)
k.showMesh()
#k.addRegion(np.array([[3,-0.5],[3.5,-0.5],[3.5,-1],[3,-1],[3,-0.5]]), 20, blocky=True, fixed=True)
#k.addRegion(np.array([[4,0],[5,0],[5,-0.5],[4,-0.5],[4,0]]), 30, blocky=True, fixed=False)

#%%
k.createSequence(params=[('wenner_alpha', [1,2,3])])
#k.createSequence()

#%%
k.forward(iplot=True, noise=0.05)
k.invert(iplot=True)

#%% show results 
# the forward initial model
k.showResults(index=0, attr='Resistivity(Ohm-m)', sens=False) # not for cR2
#k.showResults(index=0, attr='Phase(mrad)')
#k.showResults(index=0, attr='Magnitude(Ohm-m)')

# the inverted
k.showResults(index=1, attr='Resistivity(Ohm-m)', sens=False) # not for cR2

