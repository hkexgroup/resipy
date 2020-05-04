#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 11:34:56 2019

@author: jkl

"""

import numpy as np
import os
import pandas as pd
import time
import matplotlib.pyplot as plt
import resipy.meshTools as mt
from resipy.Survey import Survey
from resipy.R2 import R2

tstart = time.time()
timings = {}


testdir = 'examples/'


print('======================= GENERAL METHOD TESTS =====================')
#%% testing all importing features
k = R2()
k.createSurvey(testdir + 'dc-2d/syscal.csv', ftype='Syscal')
k.createSurvey(testdir + 'ip-2d/syscal.csv', ftype='Syscal')
k.createSurvey(testdir + 'dc-2d/syscal-normal-only.csv', ftype='Syscal')
k.addData(fname=testdir + 'dc-2d/syscal-reciprocal-only.csv', ftype='Syscal')
k.createSurvey(testdir + 'dc-2d/protocol.dat', ftype='ProtocolDC')
k.createSurvey(testdir + 'dc-3d/protocol.dat', ftype='ProtocolDC')
k.createSurvey(testdir + 'ip-2d/protocol.dat', ftype='ProtocolIP')
k.createSurvey(testdir + 'ip-3d/protocol.dat', ftype='ProtocolIP')
k.createSurvey(testdir + 'ip-3d/protocol2.dat', ftype='ProtocolIP')
k.createSurvey(testdir + 'parser/res2dinv-dd.dat', ftype='Res2Dinv')
k.createSurvey(testdir + 'parser/res2dinv-ga.dat', ftype='Res2Dinv')
k.createSurvey(testdir + 'parser/bgs-prime.dat', ftype='BGS Prime')
k.createSurvey(testdir + 'parser/sting.stg', ftype='Sting')
k.createSurvey(testdir + 'parser/abem-lund-norm.ohm', ftype='ABEM-Lund')
k.createSurvey(testdir + 'parser/Lippmann_1.tx0', ftype='Lippmann')
k.createSurvey(testdir + 'parser/Lippmann_2.tx0', ftype='Lippmann')
k.createSurvey(testdir + 'parser/Lippmann_3.tx0', ftype='Lippmann')
k.createSurvey(testdir + 'parser/syscal-new-format.csv', ftype='Syscal')

# electrode import
k = R2()
k.createSurvey(testdir + 'dc-2d-topo/syscal.csv', ftype='Syscal')
k.importElec(testdir + 'dc-2d-topo/elec.csv')

# remote detection
k = R2()
k.createSurvey(testdir + 'dc-2d-pole-dipole/syscal.csv', ftype='Syscal')
k.showPseudo(vmax=20)

# 3D survey from 2D parallel lines
k = R2(typ='R3t')
k.create3DSurvey(testdir + 'dc-2d-timelapse/data', lineSpacing=2,
                 zigzag=False, name='mergedSurvey', ftype='Syscal')

timings['methods-importing'] = time.time() - tstart


#%% filtering
k = R2()
k.createSurvey(testdir + 'ip-2d/syscal.csv')
k.filterDummy()
k.filterUnpaired() # will remove dummy but can remove more as well
k.filterElec([2])
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

k = R2()
k.createBatchSurvey(testdir + 'ip-2d-timelapse-syscal/')
k.filterRecipIP(index=-1)
k.filterRecipIP(index=-2)

timings['methods-filtering'] = time.time() - tstart


#%% error modelling
k = R2(typ='cR2')
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
# 2D flat
k = R2()
k.createSurvey(testdir + 'dc-2d/syscal.csv')
k.createMesh('quad', surface=np.array([[0, 0, 1], [3, 0, 1]]))
k.createMesh('trian', refine=1)

# 2D topo
k = R2()
k.createSurvey(testdir + 'dc-2d-topo/syscal.csv')
k.createMesh('quad')
k.createMesh('trian')

# 2D borehole (see example)

# 3D flat
#k = R2(typ='R3t') # tested in cases
#k.createSurvey(testdir + 'dc-3d/protocol.dat', ftype='ProtocolDC')
#k.createMesh()

# 3D topo
#k = R2(typ='R3t')
#k.createSurvey(testdir + 'dc-3d/protocol.dat', ftype='ProtocolDC')
#k.createMesh()

timings['methods-meshing'] = time.time() - tstart


#%% display (will be tested in the cases)
# k = R2()
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

# k = R2() # IP specific
# k.createSurvey(testdir + 'ip-2d/syscal.csv')
# k.showPseudoIP()
# #k.showHeatmap() # UI only ?
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
k = R2('resipy/invdir/test2d/', typ='R2')
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
k.createMesh(typ='trian', cl=0.1, cl_factor=5)
k.showMesh()

#k.fitErrorLin()
#k.fitErrorLME(iplot=True)
k.fitErrorPwl()
k.saveErrorData(os.path.join(k.dirname, 'dferrors.csv'))
k.saveFilteredData(os.path.join(k.dirname, 'dfdata1'), k.elec, savetyp='Res2DInv (*.dat)')
k.saveFilteredData(os.path.join(k.dirname, 'dfdata2'), k.elec, savetyp='Comma Separated Values (*.csv)')
k.err = True

k.write2in()
k.invert(modErr=True, modelDOI=True)
k.showResults(attr='Conductivity(mS/m)', doiSens=True)
k.showResults(doi=True)

# save and load project
k.saveProject(testdir + 'project')
k.loadProject(testdir + 'project.resipy')

print('elapsed: {:.4}s'.format(time.time() - t0)) # 22.7s
timings['dc-2d-topo'] = time.time() - t0


#%% test for borehole
plt.close('all')
print('-------------Testing borehole------------')
t0 = time.time()
k = R2()
k.createSurvey(testdir + 'dc-2d-borehole/protocol.dat', ftype='ProtocolDC')
#df = pd.read_csv(testdir + 'dc-2d-borehole/elec.csv')
#k.setElec(df.values[:,:3])
#buried = df.values[:,-1].astype(bool)
k.setBorehole(True)
k.showPseudo()
k.importElec(testdir + 'dc-2d-borehole/elec.csv')
k.createMesh('trian', cl=0.5, cl_factor=10)
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
k = R2(typ='cR2')
k.createSurvey(testdir + 'ip-2d/syscal.csv', ftype='Syscal')
k.showHeatmap()
k.showErrorIP()
k.dca()
k.filterManual()
k = R2(typ='cR2')
k.createSurvey(testdir + 'ip-2d/protocol.dat', ftype='ProtocolIP')
k.showPseudoIP()
k.err=True # there is already error inside the protocol.dat imported
k.invert()
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
k = R2()
k.createTimeLapseSurvey([testdir + 'dc-2d-timelapse/data/17031501.csv',
                         testdir + 'dc-2d-timelapse/data/17051601.csv',
                         testdir + 'dc-2d-timelapse/data/17040301.csv'])
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
k.fitErrorPwl()
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
print('elapsed: {:.4}s'.format(time.time() - t0))
timings['dc-2d-timelapse'] = time.time() - t0


#%% test for batch inversion with moving electrodes
plt.close('all')
print('-------------Testing Batch Inversion ------------')
t0 = time.time()
k = R2()
k.createTimeLapseSurvey(testdir + 'dc-2d-timelapse/data')
for s in k.surveys:
    print(s)
k.param['reg_mode'] = 1 # background constrained
k.param['num_xz_poly'] = 0 # need full mesh for R2.computeDiff()
for s in k.surveys:
    s.elec[3,0] = np.random.normal(s.elec[3,0], s.elec[3,0]*0.05)
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
print('elapsed: {:.4}s'.format(time.time() - t0))
timings['dc-2d-batch'] = time.time() - t0



#%% forward modelling
plt.close('all')
print('-------------Testing Forward DC Modelling ------------')
t0 = time.time()
k = R2(typ='R2')
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
k.saveSequence('resipy/invdir/seq.csv')
k.importSequence('resipy/invdir/seq.csv')
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
k = R2(typ='cR2')
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
#k = R2()
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
k = R2(typ='R3t')
k.createSurvey(testdir + 'dc-3d/protocol.dat', ftype='ProtocolDC')
k.importElec(testdir + 'dc-3d/elec.csv')
k.typ = 'R2'
k.elec2distance()

k = R2(typ='R3t')
k.createSurvey(testdir + 'dc-3d/protocol.dat', ftype='ProtocolDC')
k.importElec(testdir + 'dc-3d/elec.csv')
k.createMesh(cl=1.5, interp_method='bilinear')#, cl_factor=20, cln_factor=500)
#k.mesh.write_vtk('resipy/test/mesh3D.vtk',title='3D mesh with flat surface')
#k.err = True
k.invert(modErr=True)
k.showResults()
k.showSlice(axis='z')
k.showSlice(axis='x')
k.showSlice(axis='y')
k.showPseudoInvError()
k.showInvError()
k.saveMeshVtk()
print('elapsed: {:.4}s'.format(time.time() - t0))
timings['dc-3d'] = time.time() - t0



#%% 3D testing importing a mesh
plt.close('all')
print('-------------Testing 3D inversion ------------')
t0 = time.time()
k = R2(typ='R3t')
k.createSurvey(testdir + 'dc-3d/protocol.dat', ftype='ProtocolDC')
k.importElec(testdir + 'dc-3d/elec.csv')
k.importMesh(testdir + 'mesh/mesh3D.vtk')
k.mesh.refine() # test refining mesh 
k.addFlatError()
k.invert(modErr=False)
k.showResults() 
k.showSlice(axis='z')
k.showSlice(axis='x')
k.showSlice(axis='y')
#mesh calculations 
smesh = k.mesh.extractSurface() # this test mesh.computeNiegh as well 
tmesh = k.meshResults[0].threshold(attr='Resistivity',vmin=20,vmax=100)
print('elapsed: {:.4}s'.format(time.time() - t0))
timings['dc-3d-import-mesh'] = time.time() - t0


#%% 3D ip testing
plt.close('all')
print('-------------Testing 3D IP inversion ------------')
t0 = time.time()
k = R2(typ='cR3t')
k.createSurvey(testdir + 'ip-3d/protocol2.dat', ftype='ProtocolIP')
k.importElec(testdir + 'ip-3d/elec2.csv')
k.param['min_error'] = 0.0
k.createMesh(cl=5)
k.showMesh()
k.invert()
k.showResults()
k.showSlice(index=0)
k.showSlice(axis='z')
k.showSlice(axis='x')
k.showSlice(axis='y')
print('elapsed: {:.4}s'.format(time.time() - t0))
timings['ip-3d'] = time.time() - t0

#%% 3D column mesh -- forward modelling for 3D too 
print('----------- Testing 3D Column inversion -----------')
t0 = time.time()

k = R2(typ='R3t') # create R2 class
k.importElec(testdir + 'dc-3d-column/elec.csv') # import electrodes 
k.createMesh(typ='prism',cl=0.1,elemz=2)
a = np.array(k.mesh.elm_centre).T
idx = (a[:,1]<0.45) & (a[:,1]>-0.45) & (a[:,0]<0.45) & (a[:,0]>-0.45) # set a zone of different resistivity 
res0 = np.array(k.mesh.attr_cache['res0'])
res0[idx] = 50
k.setRefModel(res0) # set parameters for forward model 
k.showMesh(attr='res0',color_map='jet')
k.createSequence(params=[('dpdp1',1,1),('dpdp1',2,1)]) # create a sequence 
k.forward() # do forward model 

k.setRefModel(np.ones_like(res0)*100) # reset reference model 
k.invert() # invert the problem 
k.showResults(index=1) #show result 

print('elapsed: {:.4}s'.format(time.time() - t0))
timings['dc-3d-column-mesh'] = time.time() - t0

#%% print final summary information 
for key in timings.keys():
    print('{:s} : {:.2f}s'.format(key, timings[key]))
print('total time running the test = {:.4f}s'.format(time.time() - tstart))
# 

#%% test timelapse 3D -- takes a long time
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

