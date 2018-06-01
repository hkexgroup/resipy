# -*- coding: utf-8 -*-
"""
Created on Thu May 31 17:01:02 2018
Example of how to use mesh object
@author: jamyd91
"""
#import conda libraries 
import numpy as np 
#import R API libraries
import meshTools as mt
import GmshWrap as gw

#%% set up mock up geometry of an ERT survey on a slope 
surf_x=[-40, 0 , 170, 210]
surf_y=[50,50,100,100]

#%% lets pretend to put electrodes on the slope 
elec_x=np.linspace(0,170,25)#electrode x positions 
#interpolate electrode y positions 
elec_y = np.interp(elec_x,surf_x,surf_y)

#%% generate mesh of the slope 
mesh_dict = gw.tri_mesh(surf_x,surf_y,elec_x,elec_y)
#convert mesh dictionary into a mesh object
mesh=mt.mesh_obj.mesh_dict2obj(mesh_dict)

#%% show something about the mesh 
mesh.summary()
mesh.show()
#when showing the mesh we haven't actually inverted anything, so it just shows 
#the material assigned to the mesh elements by gmsh. 

#%% 
# =============================================================================
# mesh_dict2= mt.vtk_import()
# qmesh = mt.mesh_obj.mesh_dict2obj(mesh_dict2)
# 
# qmesh.summary()
# qmesh.show()
# =============================================================================

