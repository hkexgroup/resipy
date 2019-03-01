#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 15:08:19 2019

@author: jkl
"""

import numpy as np
import matplotlib.pyplot as plt
from api.R2 import R2
plt.ioff() # this disable interactive plotting

figdir = './image/paper/'

#%% 2D topo at Lancaster Castle Hill (Lancaster UK)
k = R2() # initiate an R2 instance
k.createSurvey('./api/test/syscalFileTopo.csv', ftype='Syscal') # import data
k.setElec(np.genfromtxt('./api/test/elecTopo.csv', delimiter=','))
k.pwlfit() # fit a power law
k.createMesh(typ='trian') # create quadrilateral mesh
k.invert() # run the inversion
k.showResults()

# graph section
fig, ax = plt.subplots(figsize=(7, 2))
k.showResults(ax=ax, zlim=[28, 29.6]) # show the inverted section
fig.tight_layout()
fig.savefig(figdir + 'castle.eps')
fig.show()


#%% IP Calcite precipitation inversion (Somewhere USA)
k = R2(typ = 'cR2') # initiate an R2 instance (considering there is IP data in the input data)
k.createSurvey('./api/test/IP/IP_MICP_all.csv', ftype='Syscal') # import data
k.removenested() # removing nested measurements
k.iprangefilt(0,25) # setting phase shift range to 0 < -𝟇 < 25
k.pwlfit() # adding resistance power-law error model to data
k.plotIPFit() # adding phase power-law error model to data
k.errTyp = 'pwl' # setting resistance error type flag to power-law
k.errTypIP = 'pwlip' # setting phase error type flag to power-law
k.param['a_wgt'] = 0 # "a_wgt" = 0 when there is individual resistance error
k.param['b_wgt'] = 0 # "b_wgt" = 0 when there is individual phase error
k.createMesh(typ='trian') # create triangular mesh
k.invert() # run the inversion (and write cR2.in and protocol.dat automatically)
k.showResults(attr='Sigma_real(log10)') # show the inverted real conductivity section
k.showResults(attr='Phase(mrad)') # show the inverted phase shift section

# graph section
fig, ax = plt.subplots(figsize=(6, 2))
ax.set_title('(a)')
k.showResults(attr='Sigma_real(log10)', zlim=[-8, 0], ax=ax)
fig.tight_layout()
fig.savefig(figdir + 'micp-sigma.eps')
fig.show()
fig, ax = plt.subplots(figsize=(6, 2))
ax.set_title('(b)')
k.showResults(attr='Phase(mrad)', zlim=[-8, 0], ax=ax)
fig.tight_layout()
fig.savefig(figdir + 'micp-phase.eps')
fig.show()


#%% Time-lapse RWU at Woburn (UK)
k = R2() # initiate an R2 instance
k.createTimeLapseSurvey('./api/test/testTimelapse', ftype='Syscal') # import directory with the data
k.pwlfit() # fit a power-law
k.createMesh(typ='trian', cl=0.02, cl_factor=20) # create a triangular mesh with a characteristic length of 0.5
k.showMesh()
k.invert(parallel=True) # run the inversion (and write R2.in and protocol.dat automatically)
k.showResults(index=0) # show the first inverted section
k.showResults(index=1) # show the second inverted section
k.showResults(index=1, attr='difference(percent)') # show the differences between the first and second survey

#%% graph
fig, ax = plt.subplots(figsize=(5, 2))
ax.set_title('(a) 15th March 2017')
k.showResults(ax=ax, index=1, attr='difference(percent)', vmin=0, vmax=50)
fig.tight_layout()
fig.savefig(figdir + 'woburnMarch.eps')
fig.show()

fig, ax = plt.subplots(figsize=(5, 2))
ax.set_title('(b) 27th April 2017')
k.showResults(ax=ax, index=2, attr='difference(percent)', vmin=0, vmax=50)
fig.tight_layout()
fig.savefig(figdir + 'woburnApril.eps')
fig.show()

fig, ax = plt.subplots(figsize=(5, 2))
ax.set_title('(c) 16th Mai 2017')
k.showResults(ax=ax, index=3, attr='difference(percent)', vmin=0, vmax=50)
fig.tight_layout()
fig.savefig(figdir + 'woburnMai.eps')
fig.show()


#%% ERT in River at Boxford (UK) with fixed region
k = R2()
k.createSurvey('./api/test/primeFile.dat', ftype='BGS Prime')
# following lines will add electrode position, surface points and specify if electrodes are buried or not. Similar steps are done in the GUI in (a), (b), (c)
x = np.genfromtxt('./api/test/primePosBuried.csv', delimiter=',')
k.elec[:,[0,2]] = x[:,:2] # electrode positions
surface = np.array([[0.7, 92.30],[10.3, 92.30]]) # additional surface point for the river level
buried = x[:,2].astype(bool) # specify which electrodes are buried (in the river here)
k.filterElec([21, 2]) # filter out problematic electrodes 21 and 2
k.createMesh(typ='trian', buried=buried, surface=surface, cl=0.2, cl_factor=10)
xy = k.elec[1:21,[0,2]] # adding river water level using 2 topo points
k.addRegion(xy, res0=32, blocky=True, fixed=True) # fixed river resistivity to 32 Ohm.m
k.param['b_wgt'] = 0.05 # setting up higher noise level
k.invert()
k.showResults(sens=False, vmin=1.2, vmax=2.2, zlim=[88, 93])

# graph
fig, ax = plt.subplots(figsize=(7, 2))
k.showResults(ax=ax, sens=False, vmin=1.2, vmax=2.2, zlim=[88, 93])
fig.tight_layout()
fig.savefig(figdir + 'fixedRiver.eps')
fig.show()


#%% Forward modelling for different array type



