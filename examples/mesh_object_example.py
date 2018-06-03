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
import gmshWrap as gw
import matplotlib.pyplot as plt

plt.close('all')

#%% set up mock up geometry of an ERT survey on a slope 
surf_x=[-40, 0 , 170, 210]
surf_y=[50,50,100,100]

#%% lets pretend to put electrodes on the slope 
elec_x=np.linspace(0,170,25)#electrode x positions 
#interpolate electrode y positions 
elec_y = np.interp(elec_x,surf_x,surf_y)

#axis limits 
xlim=(min(surf_x),max(surf_x))
ylim=(40,110)

#%% generate a triagnular mesh of the slope 
tri_mesh = gw.tri_mesh(surf_x,surf_y,elec_x,elec_y)

# show something about the mesh 
tri_mesh.summary()
tri_mesh.show(xlim=xlim,ylim=ylim)
#when showing the mesh we haven't actually inverted anything, so it just shows 
#the material assigned to the mesh elements by gmsh. 

#%% generate a quad mesh 
qmesh,meshx,meshy,topo,e_nodes = mt.quad_mesh(elec_x,elec_y)#
#mesh object, x node coordinates, y node coordinates, topography, x coordinate nodes with electrodes

#show the quad mesh 
qmesh.summary()
qmesh.show(xlim=xlim,ylim=ylim)

#again we havent inverted anything so the mesh is just blank 

#%% importing a vtk file 

mesh_dict = mt.vtk_import()

mesh_obj = mt.Mesh_obj.mesh_dict2obj(mesh_dict)

#show info about mesh 
mesh_obj.summary()
mesh_obj.show()
