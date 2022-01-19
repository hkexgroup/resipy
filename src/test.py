#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 11:34:56 2019

@author: jkl

"""

import numpy as np
import os
import shutil
import pandas as pd
import time
import matplotlib.pyplot as plt
import resipy.meshTools as mt
from resipy import Project, Survey
from resipy.Project import apiPath
# import warnings # too much warnings for CI
# warnings.simplefilter("always")

use_pyvista = False

tstart = time.time()
timings = {}


testdir = 'examples/'


print('======================= GENERAL METHOD TESTS =====================')
#%% testing all importing features
k = Project(typ='R2')
k.createSurvey(testdir + 'dc-2d/syscal.csv', ftype='Syscal')
k.createSurvey(testdir + 'ip-2d/syscal.csv', ftype='Syscal')
k.createSurvey(testdir + 'dc-2d/syscal-normal-only.csv', ftype='Syscal')
k.addData(fname=testdir + 'dc-2d/syscal-reciprocal-only.csv', ftype='Syscal')
k.createSurvey(testdir + 'dc-2d/protocol.dat', ftype='ProtocolDC')
k.createSurvey(testdir + 'dc-2d-borehole/xbh_prosys3.csv', ftype='Syscal')
k.createSurvey(testdir + 'dc-3d/protocol.dat', ftype='ProtocolDC')
k.createSurvey(testdir + 'ip-2d/protocol.dat', ftype='ProtocolIP')
k.createSurvey(testdir + 'ip-3d/protocol.dat', ftype='ProtocolIP')
k.createSurvey(testdir + 'ip-3d/protocol2.dat', ftype='ProtocolIP')
k.createSurvey(testdir + 'parser/res2dinv-dd.dat', ftype='ResInv')
k.createSurvey(testdir + 'parser/res2dinv-ga.dat', ftype='ResInv')
k.createSurvey(testdir + 'parser/res2dinv-multigradient.dat', ftype='ResInv')
k.createSurvey(testdir + 'parser/res2dinv-wenner32.dat', ftype='ResInv')
k.createSurvey(testdir + 'parser/bgs-prime.dat', ftype='BGS Prime')
k.createSurvey(testdir + 'parser/sting_2D_noIP.stg', ftype='Sting')
k.createSurvey(testdir + 'parser/sting_3D_noIP.stg', ftype='Sting')
k.createSurvey(testdir + 'parser/abem-lund-norm.ohm', ftype='ABEM-Lund')
k.createSurvey(testdir + 'parser/Lippmann_1.tx0', ftype='Lippmann')
k.createSurvey(testdir + 'parser/Lippmann_2.tx0', ftype='Lippmann')
k.createSurvey(testdir + 'parser/Lippmann_3.tx0', ftype='Lippmann') 
k.createSurvey(testdir + 'parser/syscal-new-format.csv', ftype='Syscal') 
k.createSurvey(testdir + 'parser/syscal_ProsysIII_IP.csv', ftype='Syscal')
k.createSurvey(testdir + 'parser/BERT_2D_topo.ohm', ftype='BERT')
k.createSurvey(testdir + 'parser/BERT_IP_2D.dat', ftype='BERT')
k.createSurvey(testdir + 'parser/DAS-1_2D_DC.data', ftype='DAS-1')
k.createSurvey(testdir + 'parser/DAS-1_3D_IPDC.data', ftype='DAS-1')
k.createSurvey(testdir + 'parser/protocolForward/R2_forward.dat', ftype='ProtocolDC')
k.createSurvey(testdir + 'parser/protocolForward/cR2_forward.dat', ftype='ProtocolIP')
k = Project(typ='R3t')
k.createSurvey(testdir + 'parser/protocolForward/R3t_forward.dat', ftype='ProtocolDC')
k.createSurvey(testdir + 'parser/protocolForward/cR3t_forward.dat', ftype='ProtocolIP')

# electrode import
k = Project(typ='R2')
k.createSurvey(testdir + 'dc-2d-topo/syscal.csv', ftype='Syscal')
k.importElec(testdir + 'dc-2d-topo/elec.csv')
k.mergeElec()

# remote detection
k = Project(typ='R2')
k.createSurvey(testdir + 'dc-2d-pole-dipole/syscal.csv', ftype='Syscal')
k.showPseudo(vmax=20)

# 3D survey from 2D parallel lines
k = Project(typ='R3t')
k.create3DSurvey(testdir + 'dc-2d-timelapse/data', lineSpacing=2,
                 zigzag=False, name='mergedSurvey', ftype='Syscal')

# 3D survey from 2D perpendicular line with one common elec
k = Project(typ='R3t')
k.create3DSurvey(testdir + 'dc-2d-pseudo3d-synthetic/data', lineSpacing=1,
                 zigzag=False, name='mergedSurvey', ftype='ProtocolDC')
k.importElec(testdir + 'dc-2d-pseudo3d-synthetic/lines-elec.csv')

timings['methods-importing'] = time.time() - tstart


#%% filtering
k = Project(typ='R2')
k.createSurvey(testdir + 'ip-2d/syscal.csv')
k.filterDummy()
k.filterUnpaired() # will remove dummy but can remove more as well
k.filterElec(['2'])
k.filterNested()
#k.filterDCA() # tested in cases
k.filterManual()
k.filterRangeIP(phimin=-10, phimax=10)
k.filterRecip(percent=20)
k.filterStack(percent=2)
k.filterRecipIP()
k.addFilteredIP()
k.filterAppResist(vmin=0, vmax=50)
k.filterTransferRes(vmin=-20, vmax=20)
k.filterZeroMeasSurveys() # remove surveys with 0 measurements (not sure it's used)
k.filterNegative() # all Tx are negative

k = Project(typ='R2')
k.createBatchSurvey(testdir + 'ip-2d-timelapse-syscal/')
k.filterRecipIP(index=-1)
k.filterRecipIP(index=-2)

timings['methods-filtering'] = time.time() - tstart


#%% error modelling
k = Project(typ='cR2')
k.createBatchSurvey(testdir + 'ip-2d-timelapse-syscal/')

k.err = True
k.write2protocol() # triggers default combined error model for DC and IP

k.showErrorIP(index=0)
k.showErrorIP(index=-2)

fig, axs = plt.subplots(5, 1, figsize=(6,6))
k.fitErrorLin(index=-1, ax=axs[0])
k.fitErrorLin(index=-2, ax=axs[1])
k.fitErrorLin(index=0, ax=axs[2])
k.fitErrorLin(index=1, ax=axs[3])
k.fitErrorLin(index=2, ax=axs[4])
fig, ax = plt.subplots()
k.fitErrorLin(index=-1, ax=ax)

k.fitErrorPwl(index=-1)
k.fitErrorPwl(index=-2)
k.fitErrorPwl(index=0)
k.fitErrorPwl(index=1)
k.fitErrorPwl(index=2)
fig, ax = plt.subplots()
k.fitErrorPwl(index=-1, ax=ax)

k.fitErrorPwlIP(index=-1)
k.fitErrorPwlIP(index=-2)
k.fitErrorPwlIP(index=0)
k.fitErrorPwlIP(index=1)
k.fitErrorPwlIP(index=2)
fig, ax = plt.subplots()
k.fitErrorPwlIP(index=-1, ax=ax)

k.fitErrorParabolaIP(index=-1)
k.fitErrorParabolaIP(index=-2)
k.fitErrorParabolaIP(index=0)
k.fitErrorParabolaIP(index=1)
k.fitErrorParabolaIP(index=2)
fig, ax = plt.subplots()
k.fitErrorParabolaIP(index=-1, ax=ax)

#k.fitErrorLME() # only tested with an R kernel

timings['methods-error-modelling'] = time.time() - tstart


#%% mesh generation (will be tested in the cases)
plt.close('all')
# 2D flat
k = Project(typ='R2')
k.createSurvey(testdir + 'dc-2d/syscal.csv')
#k.createMesh('quad', surface=np.array([[0, 0, 1], [3, 0, 1]]))
k.createMesh('trian')
k.mesh.computeNeigh()
rmesh = k.mesh.refine()
#k.mesh.connection = k.mesh.connection.astype(np.int_)
#idx = np.argsort(t) 

# 2D topo
k = Project(typ='R2')
k.createSurvey(testdir + 'dc-2d-topo/syscal.csv')
k.importElec(testdir + 'dc-2d-topo/elec.csv')
k.createMesh('quad')
k.showMesh()
k.createMesh('trian')
k.showMesh()
k.createMesh('trian',refine=1)
k.showMesh()

# 2D borehole (see example)

# 3D flat
#k = Project(typ='R3t') # tested in cases
#k.createSurvey(testdir + 'dc-3d/protocol.dat', ftype='ProtocolDC')
#k.createMesh()

# 3D topo
#k = Project(typ='R3t')
#k.createSurvey(testdir + 'dc-3d/protocol.dat', ftype='ProtocolDC')
#k.createMesh()

# 3D cylinder
# radius = 6.5/2 # cm
# angles = np.linspace(0, 2*np.pi, 13)[:-1] # radian
# celec = np.c_[radius*np.cos(angles), radius*np.sin(angles)]
# elec = np.c_[np.tile(celec.T, 8).T, np.repeat(6.5+np.arange(0, 8*5.55, 5.55)[::-1], 12)]
k = Project(typ='R3t')
# k.setElec(elec)
k.importElec(testdir + 'dc-3d-cylinder/elec.csv')
k.createMesh('cylinder', zlim=[0, 47.5], cl=0.8)
# k.importSequence(testdir + 'dc-3d-cylinder/sequence.csv')
k.createSequence([('custSeq', testdir + 'dc-3d-cylinder/sequence.csv')])
k.forward()
k.saveMesh(os.path.join(k.dirname, 'mesh.vtk'))
k.saveMesh(os.path.join(k.dirname, 'mesh.node'))
k.saveMesh(os.path.join(k.dirname, 'mesh.dat'))

# 3D tank
elec = np.array([[0,2,2],[0,2,6],[0,3,2],[0,3,6],
                  [10,2,2],[10,2,6],[10,3,2],[10,3,6],
                  [3,0,2],[5,0,2],[7,0,2],[3,0,6],[5,0,6],[7,0,6],
                  [3,5,2],[5,5,2],[7,5,2],[3,5,6],[5,5,6],[7,5,6]
                  ])
k = Project(typ='R3t')
k.setElec(elec)
k.createMesh('tank', origin=[0,0,0], dimension=[10,5,7])


# specific mesh import
# mesh = mt.tetgen_import(os.path.join(k.dirname, 'mesh.1.node'))


timings['methods-meshing'] = time.time() - tstart


#%% display (will be tested in the cases)
# k = Project(typ='R2')
# k.createTimeLapseSurvey(testdir + 'dc-2d-timelapse/data')
# k.invert(parallel=True)
# k.showPseudo(0)
# k.showPseudo(2)
# k.showError(0)
# k.showError(2)
# k.showErrorDist(0)
# k.showErrorDist(2)
# k.showInvError()
# k.showIter()
# k.showMesh()
# k.showParam()
# k.showPseudoInvError()

# k.showResults()
#k.showSlice() # in 3D cases

# k = Project(typ='R2') # IP specific
# k.createSurvey(testdir + 'ip-2d/syscal.csv')
# k.showPseudoIP()
# #k.showHeatmap() # must have k.surveys[0].filt_typ = 'Raw' or 'Filtered # TODO fix this @Sina
# k.showErrorIP()
#k.showSection() #TODO in cases or deprecate
#k.showPseudoInvErrorIP() # tested in cases

# timings['methods-showing'] = time.time() - tstart


print('elapsed for general methods test: {:.4}s'.format(time.time() - tstart))
# 77s on macbookpro corei5
print(timings)
print('=================================== cases test =======================')


#%% testing the R2 class
plt.close('all')
print('-------------Testing simple 2D inversion ------------')
t0 = time.time()
k = Project(apiPath + '/invdir/test2d/', typ='R2')
k.createSurvey(testdir + 'dc-2d-topo/syscal.csv', ftype='Syscal')
k.setTitle('Test 2D')
k.importElec(testdir + 'dc-2d-topo/elec.csv')
k.fitErrorPwl()
k.filterManual(attr='resError')
k.estimateError()
k.filterManual(attr='resError')
k.showPseudo(contour=True)
k.createMesh(typ='quad', elemx=4)
k.showMesh()
xz = [[0,2,4,6],[29.20,29.0,28.5,27.75]]
geom_input={'boundary1':xz}
k.createMesh(typ='trian', cl=0.1, cl_factor=5,geom_input=geom_input)
# mesh = mt.triMesh(k.elec['x'].values,k.elec['z'].values,geom_input=geom_input)
k.showMesh()

#k.fitErrorLin()
#k.fitErrorLME(iplot=True)
k.fitErrorPwl()
k.saveErrorData(os.path.join(k.dirname, 'dferrors.csv'))
k.saveFilteredData(os.path.join(k.dirname, 'dfdata1'), savetyp='Res2DInv (*.dat)')
k.saveFilteredData(os.path.join(k.dirname, 'dfdata2'), savetyp='Comma Separated Values (*.csv)')
k.err = True

k.invert(modErr=True, modelDOI=True)
k.showResults(attr='Conductivity(mS/m)', doiSens=True)
k.showResults(doi=True)

#filter data based on inversion error
k.filterInvError(vmin=-3, vmax=3)

# save and load project
k.saveProject(testdir + 'project')
k.loadProject(testdir + 'project.resipy')

print('elapsed: {:.4}s'.format(time.time() - t0)) # 22.7s
timings['dc-2d-topo'] = time.time() - t0


#%% test for borehole
plt.close('all')
print('-------------Testing borehole------------')
t0 = time.time()
k = Project(typ='R2')
k.createSurvey(testdir + 'dc-2d-borehole/protocol.dat', ftype='ProtocolDC')
#df = pd.read_csv(testdir + 'dc-2d-borehole/elec.csv')
#k.setElec(df.values[:,:3])
#buried = df.values[:,-1].astype(bool)
k.setBorehole(True)
k.showPseudo()
k.importElec(testdir + 'dc-2d-borehole/elec.csv')
k.createMesh('trian', cl=0.5, cl_factor=10, fmd=20)#, surface=np.array([[-5, 0],[5,0]]))
k.showMesh()
k.invert()
k.showIter(index=0)
k.showIter(index=1)
k.showResults(sens=False, contour=True)
print('elapsed: {:.4}s'.format(time.time() - t0))
timings['dc-2d-borehole'] = time.time() - t0


#%% test for IP
plt.close('all')
print('-------------Testing IP ------------')
t0 = time.time()
k = Project(typ='cR2')
k.createSurvey(testdir + 'ip-2d/syscal.csv', ftype='Syscal')
k.showHeatmap()
k.showErrorIP()
k.filterDCA()
k.filterManual()
k = Project(typ='cR2')
k.createSurvey(testdir + 'ip-2d/protocol.dat', ftype='ProtocolIP')
k.showPseudoIP()
print('k')
k.err=True # there is already error inside the protocol.dat imported
k.invert()
print('l')
k.showResults(attr='Magnitude(ohm.m)', sens=False)
k.showResults(attr='Phase(mrad)', sens=False)
k.showPseudoInvError()
k.showPseudoInvErrorIP()
# k.showInvError()
print('elapsed: {:.4}s'.format(time.time() - t0))
timings['ip-2d-topo'] = time.time() - t0
# 

#%% test for timelapse inversion
plt.close('all')
print('-------------Testing Time-lapse in // ------------')
t0 = time.time()
k = Project(apiPath + '/invdir/test2d-timelapse/')
k.createTimeLapseSurvey([testdir + 'dc-2d-timelapse/data/17031501.csv',
                         testdir + 'dc-2d-timelapse/data/17040301.csv',
                         testdir + 'dc-2d-timelapse/17051601_incomplete.csv'])
k.showPseudo(0)
k.showPseudo(2)
k.showError(0)
k.showError(2)
k.showError(-2) # combined surveys
k.showErrorDist(0)
k.showErrorDist(2)
k.showErrorDist(-2) # combined surveys
k.createMesh()
k.showMesh()
k.showParam()
k.fitErrorPwl(-1)
k.err = True
k.invert(parallel=True)
df = k.getR2out()
k.showRMS()
k.showIter()
k.showResults(index=1)
k.showResults(index=2)
k.showInvError()
k.showPseudoInvError()
k.saveInvPlots(attr='difference(percent)')

k2 = Project(apiPath + '/invdir/t/')
k2.loadResults(k.dirname)
k2.showResults()
print('elapsed: {:.4}s'.format(time.time() - t0))
timings['dc-2d-timelapse'] = time.time() - t0


#%% test for batch inversion with moving electrodes
plt.close('all')
print('-------------Testing Batch Inversion ------------')
t0 = time.time()
k = Project(typ='R2')
k.createTimeLapseSurvey(testdir + 'dc-2d-timelapse/data')
for s in k.surveys:
    print(s)
k.param['reg_mode'] = 1 # background constrained
k.param['num_xz_poly'] = 0 # need full mesh for R2.computeDiff()
for s in k.surveys:
    s.elec.loc[3, 'x'] = np.random.normal(s.elec.loc[3,'x'], s.elec.loc[3,'x']*0.05)
k.filterUnpaired()
k.createMesh()
k.fitErrorPwl()
k.err = True
k.invert(parallel=True, iMoveElec=True)
k.showResults(index=0)
k.showResults(index=1)
k.showResults(index=2, attr='difference(percent)', color_map='seismic', vmax=200)
k.saveVtks()
k.saveData('td')
shutil.rmtree('td')
print('elapsed: {:.4}s'.format(time.time() - t0))
timings['dc-2d-batch'] = time.time() - t0



#%% forward modelling
plt.close('all')
print('-------------Testing Forward DC Modelling ------------')
t0 = time.time()
k = Project(typ='R2')
k.setElec(np.c_[np.linspace(0,5.75, 24), np.zeros((24, 2))])
k.designModel(fmd=3) # interactive GUI function
k.geom_input={'polygon1':[[3, 3.5, 3.5, 3],[-0.5, -0.5, -1, -1]]}
k.createModelMesh()
k.showMesh()

k.createMesh()
k.addRegion(np.array([[1,0],[2,0],[2,-0.5],[1,-0.5],[1,0]]), 10, -3)
k.forward(iplot=True)
k.resetRegions()

# full resipy function
k.addRegion(np.array([[1,0],[2,0],[2,-0.5],[1,-0.5],[1,0]]), 10, -3)
k.addRegion(np.array([[3,-0.5],[3.5,-0.5],[3.5,-1],[3,-1],[3,-0.5]]), 20, blocky=True, fixed=False)
k.addRegion(np.array([[4,0],[5,0],[5,-0.5],[4,-0.5],[4,0]]), 30, blocky=False, fixed=True)

## full GUI function
k.createModel() # manually define 3 regions (interactive GUI function)
k.setStartingRes({0:200, 1:500, 2:20, 3:30}, {1:1, 2:2, 3:1}, {1:False, 2:False, 3:True})

# creating sequence
k.createSequence([('dpdp1', 1, 8),
                  ('dpdp2', 2, 8),
                  ('wenner_alpha', 1),
                  ('wenner_beta', 2),
                  ('wenner_gamma', 3),
                  ('schlum1', 1, 10),
                  ('schlum2', 2, 10),
                  ('multigrad', 1, 10, 2)])
k.saveSequence(k.dirname + '/seq.csv')
k.importSequence(k.dirname + '/seq.csv')
k.createSequence()
    
k.forward(iplot=True, noise=5)
# k.setRefModel([50]*k.mesh.num_elms)
k.invert()

# the forward initial model
fig, axs = plt.subplots(2, 1)
k.showResults(index=0, attr='Resistivity(Ohm-m)', sens=False, ax=axs[0]) # not for cR2
k.showResults(index=1, attr='Resistivity(Ohm-m)', sens=True, ax=axs[1], vmin=10, vmax=120) # not for cR2

print('elapsed: {:.4}s'.format(time.time() - t0))
timings['dc-2d-forward'] = time.time() - t0


#%% test forward IP modelling
plt.close('all')
print('-------------Testing Forward IP Modelling ------------')
t0 = time.time()
k = Project(typ='cR2')
k.setElec(np.c_[np.linspace(0,5.75, 24), np.zeros((24, 2))])
k.createMesh(typ='trian')
#
## full resipy function
k.addRegion(np.array([[1,0],[2,0],[2,-0.5],[1,-0.5],[1,0]]), 10, -10)

k.forward(iplot=True, noise=0.0, noiseIP=0)
k.invert(iplot=True)

# the forward initial model
fig, axs = plt.subplots(2, 1)
k.showResults(index=0, attr='Phase(mrad)', ax=axs[0])
k.showResults(index=0, attr='Magnitude(ohm.m)', ax=axs[1])

# the inverted
fig, axs = plt.subplots(2, 1)
k.showResults(index=1, attr='Phase(mrad)', ax=axs[0])
k.showResults(index=1, attr='Magnitude(ohm.m)', ax=axs[1])

print('elapsed: {:.4}s'.format(time.time() - t0))
timings['ip-2d-forward'] = time.time() - t0


#%% test Paul River (too specific) ... if so then why is it here?? - Jimmy 
#plt.close('all')
#print('-------------Testing Buried Electrodes in Fixed River ------------')
#t0 = time.time()
#k = Project(typ='R2')
#k.createSurvey(testdir + 'dc-2d-river/protocol.dat', ftype='ProtocolDC')
## following lines will add electrode position, surface points and specify if electrodes are buried or not. Similar steps are done in the GUI in (a), (b), (c)
#x = np.genfromtxt(testdir + 'dc-2d-river/elec.csv', delimiter=',')
#k.setElec(x[:,:2]) # electrode positions
#surface = np.array([[0.7, 92.30],[10.3, 92.30]]) # additional surface point for the river level
#buried = x[:,2].astype(bool) # specify which electrodes are buried (in the river here)
##k.importElec(testdir + 'dc-2d-river/elec.csv') # TODO
#k.filterElec([21, 2]) # filter out problematic electrodes 21 and 2
#k.createMesh(typ='trian', buried=buried, surface=surface, cl=0.2, cl_factor=10)
#xy = k.elec[1:21,[0,2]] # adding river water level using 2 topo points
#k.addRegion(xy, res0=32, blocky=False, fixed=False) # fixed river resistivity to 32 Ohm.m
#k.param['b_wgt'] = 0.05 # setting up higher noise level
#k.invert()
#k.showResults(sens=False, vmin=1.2, vmax=2.2, zlim=[88, 93])
#print('elapsed: {:.4}s'.format(time.time() - t0))
#timings['dc-2d-river'] = time.time() - t0


#%% 3D testing
plt.close('all')
print('-------------Testing 3D inversion ------------')
t0 = time.time()
k = Project(typ='R3t')
k.createSurvey(testdir + 'dc-3d/protocol.dat', ftype='ProtocolDC')
k.importElec(testdir + 'dc-3d/elec.csv')
k.typ = 'R2'
k.elec2distance()

k = Project(typ='R3t')
k.createSurvey(testdir + 'dc-3d/protocol.dat', ftype='ProtocolDC')
k.importElec(testdir + 'dc-3d/elec.csv')
# k.showPseudo(threed=True) # only tested in pyvista setup   
k.createMesh(cl=2, refine=1)#, interp_method='bilinear', cl_factor=20, cln_factor=500)

k.createSequence()
#k.err = True
k.invert(modErr=True)
k.showResults(use_pyvista=use_pyvista)
k.showSlice(axis='z')
k.showSlice(axis='x')
k.showSlice(axis='y')
k.showPseudoInvError()
k.showInvError()
k.saveMeshVtk()
print('elapsed: {:.4}s'.format(time.time() - t0))
timings['dc-3d'] = time.time() - t0
#timeit k.mesh.orderNodes()
#print(k.mesh)

#%% 3D testing importing and exporting a mesh 
plt.close('all')
print('-------------Testing 3D inversion with custom mesh ------------')
t0 = time.time()
k = Project(typ='R3t')
k.createSurvey(testdir + 'dc-3d/protocol.dat', ftype='ProtocolDC')
k.importElec(testdir + 'dc-3d/elec.csv')
k.importMesh(testdir + 'mesh/coarse3D.vtk')
k.mesh = k.mesh.refine() # test refining mesh 
k.addFlatError()
k.invert()

k.showResults(use_pyvista=use_pyvista)
k.showSlice(axis='z')
k.showSlice(axis='x')
k.showSlice(axis='y')
#mesh calculations 
k.mesh.exportTetgenMesh(testdir + 'mesh/tetgen_test')
smesh = k.mesh.extractSurface() # this test mesh.computeNiegh as well 
tmesh = k.meshResults[0].threshold(attr='Resistivity(ohm.m)',vmin=20,vmax=100)
print('elapsed: {:.4}s'.format(time.time() - t0))
timings['dc-3d-import-mesh'] = time.time() - t0


#%% 3D ip testing
plt.close('all')
print('-------------Testing 3D IP inversion ------------')
t0 = time.time()
k = Project(typ='cR3t')
k.createSurvey(testdir + 'ip-3d/protocol2.dat', ftype='ProtocolIP')
k.importElec(testdir + 'ip-3d/elec2.csv')
k.param['min_error'] = 0.0
k.createMesh(cl=4)
k.showMesh(use_pyvista=use_pyvista)

k.invert()

k.showResults(use_pyvista=use_pyvista)
k.showSlice(index=0)
k.showSlice(axis='z')
k.showSlice(axis='x')
k.showSlice(axis='y')
print('elapsed: {:.4}s'.format(time.time() - t0))
timings['ip-3d'] = time.time() - t0


#%% 3D column mesh 
print('----------- Testing 3D Column prism with fwd and inv -----------')
t0 = time.time()

k = Project(typ='R3t') # create R2 class
k.importElec(testdir + 'dc-3d-column-prism/elec.csv') # import electrodes 
k.createMesh(typ='prism',cl=0.1,elemz=2)

#assign regions of resistivity 
a = k.mesh.elmCentre
idx = (a[:,1]<0.45) & (a[:,1]>-0.45) & (a[:,0]<0.45) & (a[:,0]>-0.45) # set a zone of different resistivity 
res0 = np.array(k.mesh.df['res0'])
res0[idx] = 50
k.setRefModel(res0) # set parameters for forward model 

k.showMesh(attr='res0',color_map='jet',use_pyvista=use_pyvista)

#create a forward modelling sequence, bit awkward at the moment because strings need to be picked individually
xs = [0,0,1,-1]
ys = [-1,1,0,0]
zs = np.unique(k.elec['z'].values)
seqIdx = []
for i in range(4): # makes down column strings
    sb = (k.elec['x'] == xs[i]) & (k.elec['y'] == ys[i]) # bool on string
    si = [] # index on string
    for j in range(len(sb)):
        if sb[j]:
            si.append(j)
    seqIdx.append(si)
for i in range(len(zs)):
    sb = (k.elec['z'] == zs[i]) # bool on string
    si = [] # index on string
    for j in range(len(sb)):
        if sb[j]:
            si.append(j)
    seqIdx.append(si)
    

k.createSequence(params=[('dpdp1',1,1),('dpdp1',2,1)], seqIdx=seqIdx) # create a sequence 
k.forward() # do forward model 

k.setRefModel(np.ones_like(res0)*100) # reset reference model 
k.invert() # invert the problem 
k.showResults(index=1, use_pyvista=use_pyvista) #show result 

print('elapsed: {:.4}s'.format(time.time() - t0))
timings['dc-3d-column-mesh'] = time.time() - t0


#%% test 3D column inversion on tetrahedral mesh
print('---------- Testing 3D columns tetrahedral ------------')
t0 = time.time()
k = Project(typ='R3t')
k.createSurvey(testdir + 'dc-3d-timelapse-column/protocol.dat', ftype='ProtocolDC')
k.importElec(testdir + 'dc-3d-timelapse-column/elec.csv')
k.importMesh(testdir + 'dc-3d-timelapse-column/mesh.msh')
k.param['num_xy_poly'] = 0
k.param['zmin'] = -np.inf
k.param['zmax'] = np.inf
k.invert()
k.showResults(use_pyvista=use_pyvista)

print('elapsed: {:.4}s'.format(time.time() - t0))
timings['dc-3d-column-mesh'] = time.time() - t0


#%% test timelapse 3D -- takes a long time
print('----------- Testing 3D time-lapse inversion -----------')
t0 = time.time()
print('------------- 3D time-lapse difference inversion ---------------')
k = Project(typ='R3t')
k.createTimeLapseSurvey(testdir + 'dc-3d-timelapse-protocol/data' ,ftype='ProtocolDC')
k.importElec(testdir + 'dc-3d-timelapse-protocol/elec/electrodes3D-1.csv')
k.createMesh()
k.invert()

k.showResults(index=0,use_pyvista=use_pyvista)
k.showResults(index=1,use_pyvista=use_pyvista)

k.mesh.orderNodes()
print('elapsed: {:.4}s'.format(time.time() - t0))
timings['dc-3d-timelapse'] = time.time() - t0

#%% test pseudo 3D inversion
print('----------- Testing pseudo 3D inversion -----------')
t0 = time.time()
k = Project(typ='R2')
k.createPseudo3DSurvey(testdir + 'dc-2d-pseudo3d-synthetic/data', lineSpacing=1,
                 ftype='ProtocolDC')

## manually setting up electrodes
## rotating middle electrodes line here
# rotmat = np.array([[np.cos(0.2), np.sin(0.2)], 
#                     [-np.sin(0.2), np.cos(0.2)]])
# xy = np.array([k.elec.loc[24:47,'x'].values, k.elec.loc[24:47,'y'].values])
# newmat = np.dot(rotmat, xy).T
# elecTemp = k.elec.copy()
# elecTemp.loc[0:23, 'y'] = 1
# elecTemp.loc[24:47,'x'] = newmat[:,0].copy()*0.6 + 2
# elecTemp.loc[24:47,'y'] = newmat[:,1].copy()*2 - 1
# elecTemp.loc[24:35,'z'] = np.linspace(0,2,12)
# elecTemp.loc[36:47,'z'] = np.linspace(2,1,12)
# k.pseudo3DSurvey.elec = elecTemp.copy()
# k._updatePseudo3DSurvey()

## or load the files with 3D-like labels for elec positions of all lines
k.importPseudo3DElec(testdir + 'dc-2d-pseudo3d-synthetic/lines-elec.csv')
k.createMultiMesh(typ='trian', runParallel=True)
# k.showPseudo3DMesh(cropMesh=True) # only works with pyvista - thus commented for test
k.invertPseudo3D(runParallel=True)
# k.showResults(index=-1, cropMesh=True, clipCorners=False, pseudo3DContour=True) # only works with pyvista - thus commented for test
print('elapsed: {:.4}s'.format(time.time() - t0))
timings['dc-2d-pseudo3d'] = time.time() - t0

#%% print final summary information 
for key in timings.keys():
    print('{:s} : {:.2f}s'.format(key, timings[key]))
print('total time running the test = {:.4f}s'.format(time.time() - tstart))

