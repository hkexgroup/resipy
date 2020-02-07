#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 11:34:56 2019

@author: jkl

RUN ALL SECTION ON THE TEST AND CHECK THE GRAPH PRODUCED

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
k.createSurvey(testdir + 'dc-2d/protocol.dat', ftype='Protocol')
k.createSurvey(testdir + 'dc-3d/protocol.dat', ftype='Protocol')
k.createSurvey(testdir + 'ip-2d/protocol.dat', ftype='ProtocolIP')
#k.createSurvey(testdir + 'ip-3d/protocol.dat', ftype='ProtocolIP')
k.createSurvey(testdir + 'parser/res2dinv-dd.dat', ftype='Res2Dinv')
k.createSurvey(testdir + 'parser/res2dinv-ga.dat', ftype='Res2Dinv')
#k.createSurvey(testdir + 'parsers/bgs-prime.csv', ftype='BGS Prime')
k.createSurvey(testdir + 'parser/sting.stg', ftype='Sting')
#k.createSurvey(testdir + 'parsers/abem-lund.csv', ftype='ABEM-Lund')
#k.createSurvey(testdir + 'parsers/lippmann.csv', ftype='Lippmann')

# electrode import
k = R2()
k.createSurvey(testdir + 'dc-2d-topo/syscal.csv')
k.importElec(testdir + 'dc-2d-topo/elec.csv')

k = R2(typ='R3t')
k.create3DSurvey(testdir + 'dc-2d-timelapse/data')

timings['methods-importing'] = time.time() - tstart

#%% filtering
k = R2()
k.createSurvey(testdir + 'ip-2d/syscal.csv')
k.filterDummy()
k.filterUnpaired() # will remove dummy but can remove more as well
k.filterElec([2])
#k.filterNegative() # all Tx are negative
k.filterNested()
#k.filterDCA() # tested in cases
k.filterManual()
k.filterRangeIP(phimin=-10, phimax=10)
k.filterRecip(percent=20)
k.filterRecipIP()
k.filterZeroMeasSurveys() # remove surveys with 0 measurements (not sure it's used)

timings['methods-filtering'] = time.time() - tstart

#%% error modelling
k = R2()
k.createBatchSurvey(testdir + 'ip-2d-timelapse-syscal/')
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

timings['methods-error-modelling'] = time.time() - tstart

#%% mesh generation (will be tested in the cases)
# 2D flat
k = R2()
k.createSurvey(testdir + 'dc-2d/syscal.csv')
k.createMesh('quad')
k.createMesh('trian')

# 2D topo
k = R2()
k.createSurvey(testdir + 'dc-2d-topo/syscal.csv')
k.createMesh('quad')
k.createMesh('trian')

# 2D borehole (see example)

# 3D flat
#k = R2(typ='R3t') # tested in cases
#k.createSurvey(testdir + 'dc-3d/protocol.dat', ftype='Protocol')
#k.createMesh()

# 3D topo
#k = R2(typ='R3t')
#k.createSurvey(testdir + 'dc-3d/protocol.dat', ftype='Protocol')
#k.createMesh()

timings['methods-meshing'] = time.time() - tstart

#%% display (will be tested in the cases)
k = R2()
k.createTimeLapseSurvey(testdir + 'dc-2d-timelapse/data')
k.invert(parallel=True)
k.showPseudo(0)
k.showPseudo(2)
k.showError(0)
k.showError(2)
k.showErrorDist(0)
k.showErrorDist(2)
k.showInvError()
k.showIter()
k.showMesh()
k.showParam()
k.showPseudoInvError()

k.showResults()
#k.showSlice() # in 3D cases

k = R2() # IP specific
k.createSurvey(testdir + 'ip-2d/syscal.csv')
k.showPseudoIP()
#k.showHeatmap() # UI only ?
k.showErrorIP()
#k.showSection() #TODO in cases or deprecate
#k.showPseudoInvErrorIP() # tested in cases

timings['methods-showing'] = time.time() - tstart

print('elapsed for general methods test: {:.4}s'.format(time.time() - tstart))
# 77s on macbookpro corei5
print(timings)
print('=================================== cases test =======================')


#%% testing the R2 class
plt.close('all')
print('-------------Testing simple 2D inversion ------------')
t0 = time.time()
k = R2()
k.createSurvey(testdir + 'dc-2d-topo/syscal.csv', ftype='Syscal')
k.showPseudo(contour=True)
k.importElec(testdir + 'dc-2d-topo/elec.csv')
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
print('elapsed: {:.4}s'.format(time.time() - t0)) # 22.7s
timings['dc-2d-topo'] = time.time() - t0

#%% test for borehole
plt.close('all')
print('-------------Testing borehole------------')
t0 = time.time()
k = R2()
k.createSurvey(testdir + 'dc-2d-borehole/protocol.dat', ftype='Protocol')
#df = pd.read_csv(testdir + 'dc-2d-borehole/elec.csv')
#k.setElec(df.values[:,:3])
#buried = df.values[:,-1].astype(bool)
k.importElec(testdir + 'dc-2d-borehole/elec.csv')
k.createMesh('trian', cl=0.5, cl_factor=20)
k.createMesh('quad', elemx=12)
k.invert()
k.showIter(index=0)
k.showIter(index=1)
k.showResults(sens=False, contour=True)
print('elapsed: {:.4}s'.format(time.time() - t0))
timings['dc-2d-borehole'] = time.time() - t0


#%% test model DOI
plt.close('all')
print('============== Test modelDOI() =================')
k = R2()
#k.createSurvey('../examples/dc-2d-topo/syscal.csv')
#k.importElec('../examples/dc-2d-topo/elec.csv')
#k.createSurvey(testdir + 'ip-2d/protocol.dat', ftype='ProtocolIP')
k.createSurvey('../examples/dc-2d/syscal.csv')
k.createMesh('trian')
k.invert(modelDOI=True)
k.showResults(doiSens=True)
k.showResults(doi=True)


#%% test for IP
plt.close('all')
print('-------------Testing IP ------------')
t0 = time.time()

k = R2(typ='cR2')
k.createSurvey(testdir + 'ip-2d/protocol.dat', ftype='ProtocolIP')
k.err=True # there is already error inside the protocol.dat imported
k.invert()
k.showResults(attr='Magnitude(Ohm-m)', sens=False)
k.showResults(attr='Phase(mrad)', sens=False)
k.showPseudoInvError()
k.showPseudoInvErrorIP()
k.showInvError()
print('elapsed: {:.4}s'.format(time.time() - t0))
timings['ip-2d-topo'] = time.time() - t0


#%% test for timelapse inversion
plt.close('all')
print('-------------Testing Time-lapse in // ------------')
t0 = time.time()
k = R2()
k.createTimeLapseSurvey([testdir + 'dc-2d-timelapse/data/17031501.csv',
                         testdir + 'dc-2d-timelapse/data/17051601.csv',
                         testdir + 'dc-2d-timelapse/data/17040301.csv'])
#k.createTimeLapseSurvey(testdir + 'dc-2d-timelapse/data') # dirname or list of files
k.fitErrorPwl()
k.err = True
k.invert(iplot=False, parallel=True)#, modelDOI=True)
k.saveInvPlots(attr='difference(percent)')
k.showResults(index=1)
k.showResults(index=2)
print('elapsed: {:.4}s'.format(time.time() - t0))
timings['dc-2d-timelapse'] = time.time() - t0


#%% test for batch inversion with moving electrodes
plt.close('all')
print('-------------Testing Batch Inversion ------------')
t0 = time.time()
k = R2()
k.createBatchSurvey(testdir + 'dc-2d-timelapse/data')
for s in k.surveys:
    s.elec[3,0] = np.random.normal(s.elec[3,0], s.elec[3,0]*0.05)
k.filterUnpaired()
k.createMesh()
k.fitErrorPwl()
k.err = True
k.invert(parallel=True, iMoveElec=True)
k.showResults(index=0)
k.showResults(index=1)
k.showResults(index=2)
print('elapsed: {:.4}s'.format(time.time() - t0))
timings['dc-2d-batch'] = time.time() - t0



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
timings['ip-2d-forward'] = time.time() - t0


#%% test Paul River (too specific)
#plt.close('all')
#print('-------------Testing Buried Electrodes in Fixed River ------------')
#t0 = time.time()
#k = R2()
#k.createSurvey(testdir + 'dc-2d-river/protocol.dat', ftype='Protocol')
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
k.createSurvey(testdir + 'dc-3d/protocol.dat', ftype='Protocol')
k.importElec(testdir + 'dc-3d/elec.csv')
k.createMesh(cl=1.5, interp_method='bilinear')#, cl_factor=20, cln_factor=500)
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
timings['dc-3d'] = time.time() - t0



#%% 3D testing importing a mesh
plt.close('all')
print('-------------Testing 3D inversion ------------')
t0 = time.time()
k = R2(typ='R3t')
k.createSurvey(testdir + 'dc-3d/protocol.dat', ftype='Protocol')
k.importElec(testdir + 'dc-3d/elec.csv')
k.importMesh(testdir + 'mesh/mesh3D.vtk')
#k.computeModelError()
#k.write2in()
#k.param = param
k.invert()
k.showResults() 
k.showSlice(axis='z')
k.showSlice(axis='x')
k.showSlice(axis='y')
print('elapsed: {:.4}s'.format(time.time() - t0))
timings['dc-3d-import-mesh'] = time.time() - t0


#%% 3D ip testing
plt.close('all')
print('-------------Testing 3D IP inversion ------------')
t0 = time.time()
k = R2(typ='cR3t')
k.createSurvey(testdir + 'ip-3d/protocol.dat', ftype='Protocol')
k.importElec(testdir + 'ip-3d/elec.csv')
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
timings['ip-3d'] = time.time() - t0


#%% 3D with moving electrodes --takes a long time
plt.close('all')
print('-------------Testing 3D inversion ------------')
t0 = time.time()
k = R2(typ='R3t')
#k.createBatchSurvey('resipy/test/3d/data/', ftype='Protocol')
k.createTimeLapseSurvey(testdir + 'dc-3d-timelapse-protocol/data/', ftype='Protocol')
elecList = [np.genfromtxt(testdir + 'dc-3d-timelapse-protocol/elec/' + f, delimiter=',') for f in os.listdir(testdir + 'dc-3d-timelapse-protocol/elec/')]
k.setElec(elec=None, elecList=elecList)
k.createMesh(cl=2)
k.param['reg_mode'] = 1 # background regularization
k.err = True # test using estimated error model 
k.errTyp = 'survey'
k.estimateError()
k.computeModelError()
k.invert(parallel=True, iMoveElec=True)
k.showInParaview()
print('elapsed: {:.4}s'.format(time.time() - t0))
timings['ip-3d-timelapse-moving-elec'] = time.time() - t0

print(timings)
print('total time running the test = {:.4f}s'.format(time.time() - tstart))


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

