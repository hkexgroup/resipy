#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 17:44:41 2018
3d mesh example 
@author: jimmy
"""
#import modules
import numpy as np
import pandas as pd 
import sys
sys.path.append('..\src')#add pyr2 location
import api.meshTools as mt

#%% 
data = pd.read_csv(r"../src/api/test/init_elec_locs.csv")#electrode position file

#reassign the lists to numpy arrays
elec_x = np.array(data['x'])
elec_y = np.array(data['y'])
elec_z = np.array(data['z'])+2
amesh = mt.tetra_mesh(elec_x,elec_y,elec_z, interp_method = 'idw', cl = 2)
amesh.write_vtk('mesh3D.vtk')
amesh.show(alpha=1,xlim=[min(elec_x),max(elec_x)],ylim=[min(elec_y),max(elec_y)],zlim=[min(elec_z),max(elec_z)])
mt.points2vtk(elec_x,elec_y,elec_z,file_name='electrodes.vtk')
