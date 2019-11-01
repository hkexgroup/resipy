#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 11:34:56 2019

@author: jkl

RUN ALL SECTION ON THE TEST AND CHECK THE GRAPH PRODUCED

"""

import numpy as np
import os
import time
import matplotlib.pyplot as plt
import resipy.meshTools as mt
from resipy.Survey import Survey
from resipy.R2 import R2

tstart = time.time()

#%% very simple example
k = R2()
k.createSurvey('resipy/test/IP/syscalFileIP.csv')
#k.createModel()
#k.invert()
#k.saveInvPlots()
#k.pseudoError()

#%% testing the R2 class
plt.close('all')
print('-------------Testing simple 2D inversion ------------')
t0 = time.time()
k = R2()
k.createSurvey('./resipy/test/syscalFileTopo.csv', ftype='Syscal')
k.showPseudo(contour=True)
k.importElec('./resipy/test/elecTopo.csv')
#k.createMesh(typ='quad',elemx=4)
#k.showMesh()
k.createMesh(typ='trian',cl=0.1, cl_factor=5)
k.showMesh()
#k.fitErrorLin()
k.fitErrorPwl()
k.err = True
#k.fitErrorLME(iplot=True)
k.computeModelError()

k.write2in()
#k.write2protocol(err=True, errTot=True)
k.invert()
k.showResults(attr='Conductivity(mS/m)')
#k.showInParaview()
print('elapsed: {:.4}s'.format(time.time() - t0))

#%% testing error models
k = R2()
k.createBatchSurvey('./resipy/test/timelapseIP')
k.fitErrorLin(index=-1)
k.fitErrorLin(index=-2)
k.fitErrorLin(index=0)
k.fitErrorLin(index=1)
k.fitErrorLin(index=2)

k.fitErrorPwl(index=-1)
k.fitErrorPwl(index=-2)
k.fitErrorPwl(index=0)
k.fitErrorPwl(index=1)
k.fitErrorPwl(index=2)

k.fitErrorPwlIP(index=-1)
k.fitErrorPwlIP(index=-2)
k.fitErrorPwlIP(index=0)
k.fitErrorPwlIP(index=1)
k.fitErrorPwlIP(index=2)

k.fitErrorParabolaIP(index=-1)
k.fitErrorParabolaIP(index=-2)
k.fitErrorParabolaIP(index=0)
k.fitErrorParabolaIP(index=1)
k.fitErrorParabolaIP(index=2)


#%% test for borehole
plt.close('all')
print('-------------Testing borehole------------')
t0 = time.time()
k = R2()
k.createSurvey('./resipy/test/protocolXbh.dat', ftype='forwardProtocolDC')
x = np.genfromtxt('./resipy/test/elecXbh.csv', delimiter=',')
k.elec[:,[0,2]] = x[:,:2]
buried = x[:,2].astype(bool)
k.createMesh('trian', buried=buried, cl=0.5, cl_factor=20)
k.showMesh()
k.createMesh('quad', buried=buried, elemx=12)
k.showMesh()
k.invert()
k.showIter(index=0)
k.showIter(index=1)
k.showResults(sens=False, contour=True)
print('elapsed: {:.4}s'.format(time.time() - t0))


#%% test for IP
plt.close('all')
print('-------------Testing IP ------------')
t0 = time.time()

k = R2(typ='cR2')
#k.createSurvey('resipy/test/IP/rifleday8.csv', ftype='Syscal')
k.createSurvey('resipy/test/IP/syscalFileIP.csv')
k.showPseudoIP()
k.showErrorIP()
k.fitErrorPwlIP()
k.fitErrorParabolaIP()

k = R2(typ='cR2')
k.createSurvey('resipy/test/IP/protocolIP2D.dat', ftype='ProtocolIP')
k.err=True # there is already error inside the protocol.dat imported
k.invert()
k.showResults(attr='Magnitude(Ohm-m)', sens=False)
k.showResults(attr='Phase(mrad)', sens=False)
k.showPseudoInvError()
k.showPseudoInvErrorIP()
k.showInvError()
print('elapsed: {:.4}s'.format(time.time() - t0))


#%% test for timelapse inversion
plt.close('all')
print('-------------Testing Time-lapse in // ------------')
t0 = time.time()
k = R2()
k.createTimeLapseSurvey(['resipy/test/testTimelapse/17031501.csv',
                         'resipy/test/testTimelapse/17051601.csv',
                         'resipy/test/testTimelapse/17040301.csv'])
#k.createTimeLapseSurvey('resipy/test/testTimelapse')
k.fitErrorLin()
k.fitErrorPwl()
k.err = True
k.invert(iplot=False, parallel=True, ncores=2)
k.saveInvPlots(attr='difference(percent)')
k.showResults(index=1)
k.showResults(index=2)
#k.showResults(index=3) # file 17051601.csv is removed from testTimelapse folder
print('elapsed: {:.4}s'.format(time.time() - t0))


#%% test for timelapse inversion (reg_mode=1)
plt.close('all')
print('-------------Testing Time-lapse in // ------------')
t0 = time.time()
k = R2()
k.createTimeLapseSurvey('resipy/test/testTimelapse')
k.fitErrorLin()
k.fitErrorPwl()
k.err = True
k.param['reg_mode'] = 1
k.invert(iplot=False, parallel=True)
k.showResults(index=1)
k.showResults(index=2)
#k.showResults(index=3) # file 17051601.csv is removed from testTimelapse folder
print('elapsed: {:.4}s'.format(time.time() - t0))


#%% test for batch inversion with moving electrodes
plt.close('all')
print('-------------Testing Batch Inversion ------------')
from resipy.R2 import R2
t0 = time.time()
k = R2()
k.createBatchSurvey('resipy/test/testTimelapse')
for s in k.surveys:
    s.elec[3,0] = np.random.normal(s.elec[3,0], s.elec[3,0]*0.05)
k.filterUnpaired()
k.filterElec([])
k.createMesh('trian')
k.fitErrorPwl()
k.err = True
k.invert(parallel=True, iMoveElec=True)
k.showResults(index=0)
k.showResults(index=1)
k.showResults(index=2)
#k.showResults(index=3) # file 17051601.csv is removed from testTimelapse folder
print('elapsed: {:.4}s'.format(time.time() - t0))


#%% test mesh with buried electrodes
#
#print('-------------Testing Buried electrodes Inversion ------------')
#k = R2()
#k.createSurvey('resipy/test/syscalFile.csv')
#elec = k.elec
#elec[:,2] = -0.5 # let's bury them all
##elec[0,2] = 0
#k.setElec(elec)
#buried = np.ones(elec.shape[0], dtype=bool)
##buried[[0,-1]] = False # comment it and it will work
#surface = np.array([[0,0],[7,0]])
#k.createMesh(typ='quad', buried=buried, surface=surface)
#k.showMesh()

#k = R2()
#k.createSurvey('resipy/test/syscalFile.csv', ftype='Syscal')
#k.elec[:,2] = np.tile(np.arange(0,-12,-1),2)
#k.elec[:,0] = np.repeat([0,8], 12)
#k.elec[1:11,:2] = np.c_[np.ones(10)*2, np.linspace(0,-4,10)]
#buried = np.ones(k.elec.shape[0], dtype=bool)
#buried[[0,12]] = False
#surface = np.array([[-2,0],[10,0]])
#k.createMesh(typ='quad')
#k.showMesh()
#k.createMesh(typ='trian', buried=buried, cl=0.5, cl_factor=5) # works well
#k.showMesh()
#k.invert()
#k.showResults()
#
##%%
#k = R2()
#k.createSurvey('./resipy/test/syscalFile.csv')
#k.elec[3,1] = -1
#buried = np.zeros(k.elec.shape[0], dtype=bool)
#buried[3] = True
#k.createMesh('quad', buried=buried)
#k.invert()

#%% forward modelling
plt.close('all')
print('-------------Testing Forward DC Modelling ------------')
t0 = time.time()
k = R2(typ='R2')
k.setElec(np.c_[np.linspace(0,5.75, 24), np.zeros((24, 2))])
k.createMesh(typ='trian')

#k.createSequence(params=[('wenner_alpha',1), # uncomment for wenner array
#                         ('wenner_alpha',2),
#                         ('wenner_alpha',3),
#                         ('wenner_alpha',4),
#                         ('wenner_alpha',5),
#                         ('wenner_alpha',6),
#                         ('wenner_alpha',7),
#                         ('wenner_alpha',8),
#                         ('wenner_alpha',9),
#                         ('wenner_alpha',10)])

## full resipy function
k.addRegion(np.array([[1,0],[2,0],[2,-0.5],[1,-0.5],[1,0]]), 10, -3)
k.addRegion(np.array([[3,-0.5],[3.5,-0.5],[3.5,-1],[3,-1],[3,-0.5]]), 20, blocky=True, fixed=True)
k.addRegion(np.array([[4,0],[5,0],[5,-0.5],[4,-0.5],[4,0]]), 30, blocky=True, fixed=False)

## full GUI function
#k.createModel() # manually define 3 regions
#k.assignRes0({1:500, 2:20, 3:30}, {1:1, 2:2, 3:1}, {1:False, 2:False, 3:True})

# creating sequence
k.createSequence([('dpdp1', 1, 8),
                  ('wenner_alpha', 1),
                  ('wenner_alpha', 2)])
    

k.forward(iplot=True, noise=0.05)
k.invert(iplot=True)

# the forward initial model
#k.showResults(index=0, attr='Resistivity(Ohm-m)', sens=False) # not for cR2
#k.showResults(index=0, attr='Phase(mrad)')
#k.showResults(index=0, attr='Magnitude(Ohm-m)')

# the inverted
k.showResults(index=1, attr='Resistivity(Ohm-m)', sens=True, vmin=10, vmax=120) # not for cR2
#k.showResults(index=1, attr='Phase(mrad)')
#k.showResults(index=1, attr='Magnitude(Ohm-m)')
print('elapsed: {:.4}s'.format(time.time() - t0))


#%% test forward IP modelling
plt.close('all')
print('-------------Testing Forward IP Modelling ------------')
t0 = time.time()
k = R2(typ='cR2')
k.elec = np.c_[np.linspace(0,5.75, 24), np.zeros((24, 2))]
k.createMesh(typ='trian')
#
## full resipy function
k.addRegion(np.array([[1,0],[2,0],[2,-0.5],[1,-0.5],[1,0]]), 10, -3)

k.forward(iplot=True, noise=0.05, noiseIP=1)
k.invert(iplot=True)

# the forward initial model
k.showResults(index=0, attr='Phase(mrad)')
k.showResults(index=0, attr='Magnitude(Ohm-m)')

# the inverted
k.showResults(index=1, attr='Phase(mrad)')
k.showResults(index=1, attr='Magnitude(Ohm-m)')
print('elapsed: {:.4}s'.format(time.time() - t0))


#%% test Paul River
plt.close('all')
from resipy.R2 import R2
print('-------------Testing Buried Electrodes in Fixed River ------------')
t0 = time.time()
k = R2()
k.createSurvey('./resipy/test/river-protocol.dat', ftype='Protocol')
# following lines will add electrode position, surface points and specify if electrodes are buried or not. Similar steps are done in the GUI in (a), (b), (c)
x = np.genfromtxt('./resipy/test/river-elec.csv', delimiter=',')
k.setElec(x[:,:2]) # electrode positions
surface = np.array([[0.7, 92.30],[10.3, 92.30]]) # additional surface point for the river level
buried = x[:,2].astype(bool) # specify which electrodes are buried (in the river here)
k.filterElec([21, 2]) # filter out problematic electrodes 21 and 2
k.createMesh(typ='trian', buried=buried, surface=surface, cl=0.2, cl_factor=10)
xy = k.elec[1:21,[0,2]] # adding river water level using 2 topo points
k.addRegion(xy, res0=32, blocky=False, fixed=False) # fixed river resistivity to 32 Ohm.m
k.param['b_wgt'] = 0.05 # setting up higher noise level
k.invert()
k.showResults(sens=False, vmin=1.2, vmax=2.2, zlim=[88, 93])
print('elapsed: {:.4}s'.format(time.time() - t0))


#%% 3D testing
from resipy.R2 import R2
plt.close('all')
print('-------------Testing 3D inversion ------------')
t0 = time.time()
k = R2(typ='R3t')
k.createSurvey('resipy/test/protocol3D.dat', ftype='Protocol')
elec = np.genfromtxt('resipy/test/electrodes3D.csv',delimiter=',')
k.setElec(elec)
k.createMesh(cl=1.5,interp_method='bilinear')#, cl_factor=20, cln_factor=500)
#k.mesh.write_vtk('resipy/test/mesh3D.vtk',title='3D mesh with flat surface')
#k.computeModelError()
#k.err = True
k.invert()
k.showResults() 
k.showSlice(axis='z')
k.showSlice(axis='x')
k.showSlice(axis='y')
k.showPseudoInvError()
k.showInvError()
print('elapsed: {:.4}s'.format(time.time() - t0))


#%% 3D testing importing a mesh
from resipy.R2 import R2
plt.close('all')
print('-------------Testing 3D inversion ------------')
t0 = time.time()
k = R2(typ='R3t')
k.createSurvey('resipy/test/protocol3D.dat', ftype='Protocol')
elec = np.genfromtxt('resipy/test/electrodes3D.csv',delimiter=',')
k.setElec(elec)
k.importMesh('resipy/test/mesh3D.vtk')
#k.computeModelError()
#k.write2in()
#k.param = param
k.invert()
k.showResults() 
k.showSlice(axis='z')
k.showSlice(axis='x')
k.showSlice(axis='y')
print('elapsed: {:.4}s'.format(time.time() - t0))


#%% 3D ip testing
plt.close('all')
print('-------------Testing 3D IP inversion ------------')
t0 = time.time()
k = R2(typ='cR3t')
k.createSurvey('resipy/test/IP/protocol3Dip.dat', ftype='Protocol')
elec = np.genfromtxt('resipy/test/electrodes3Dip.csv', delimiter=',')
k.setElec(elec)
k.createMesh(cl=3)
k.showMesh()
k.invert()
k.showResults()
k.showSlice(index=0)
k.showSlice(axis='z')
k.showSlice(axis='x')
k.showSlice(axis='y')
k.showInParaview()
print('elapsed: {:.4}s'.format(time.time() - t0))


#%% 3D with moving electrodes (specialy dedicated to Jimmy ;)
from resipy.R2 import R2
plt.close('all')
print('-------------Testing 3D inversion ------------')
t0 = time.time()
k = R2(typ='R3t')
#k.createBatchSurvey('resipy/test/3d/data/', ftype='Protocol')
k.createTimeLapseSurvey('resipy/test/3d/data/', ftype='Protocol')
elecList = [np.genfromtxt('resipy/test/3d/elec/' + f, delimiter=',') for f in os.listdir('resipy/test/3d/elec/')]
k.setElec(elec=None, elecList=elecList)
k.createMesh(cl=2)
k.param['reg_mode'] = 1 # background regularization
k.err = True # test using estimated error model 
k.errTyp = 'survey'
k.estimateError()
k.computeModelError()
k.invert(parallel=True, iMoveElec=True)
k.showInParaview()
#print('elapsed: {:.4}s'.format(time.time() - t0))



print('total time running the test = {:.4f}s'.format(time.time() - tstart))


#%% test timelapse 3D -- takes a long time
#from resipy.R2 import R2
#k = R2(typ='R3t')
#k.createBatchSurvey('resipy/test/timelapse3D/dataLeadingRidge')
#k.importElec('resipy/test/timelapse3D/elecLeadingRidge.csv')
#k.pseudo()
#k.pwlfit()
#k.createMesh('tetra', cl=0.3, cl_factor=5)
#k.showMesh()
#k.invert(parallel=True)
#k.showResults()
#k.showInParaview()

