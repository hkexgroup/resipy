# -*- coding: utf-8 -*-
"""
Created on Thu May 31 17:01:02 2018
Example of how to use mesh object
@author: jamyd91
"""
#import conda libraries 
import numpy as np
import sys 
#import R2 API libraries
sys.path.append(r"../src")
import api.meshTools as mt
import matplotlib.pyplot as plt

plt.close('all')

#%% set up mock up geometry of an ERT survey on a slope 
surf_x=[-40, 0 , 170, 210]
surf_y=[50,50,100,100]

#%% lets pretend to put electrodes on the slope 
elec_x=np.linspace(0,170,25)#electrode x positions 
#interpolate electrode y positions 
elec_z = np.interp(elec_x,surf_x,surf_y)

#axis limits 
xlim=(min(surf_x),max(surf_x))
ylim=(40,110)

#%% generate a triagnular mesh of the slope 
geom_input = {'surface':[[-40,210],[50,100]],}
tri_mesh= mt.tri_mesh(elec_x,elec_z, geom_input=geom_input)

# show something about the mesh 
tri_mesh.summary()
tri_mesh.show(xlim=xlim,ylim=ylim)
tri_mesh.write_vtk('tri_mesh.vtk')
#when showing the mesh we haven't actually inverted anything, so it just shows 
#the material assigned to the mesh elements by gmsh. 

#%% generate a quad mesh 
qmesh,meshx,meshy,topo,e_nodes = mt.quad_mesh(elec_x,elec_z)#
#mesh object, x node coordinates, y node coordinates, topography, x coordinate nodes with electrodes

#show the quad mesh 
qmesh.summary()
qmesh.show(xlim=xlim,ylim=ylim)
qmesh.write_vtk('qmesh.vtk')

#again we havent inverted anything so the mesh is just blank 

#%% importing a vtk file 

mesh_obj = mt.vtk_import('qmesh.vtk')

#show info about mesh 
mesh_obj.summary()
mesh_obj.show(xlim=xlim,zlim=ylim)

#%% mesh with boreholes 
# import gmshWrap as gw
# 
surf_x = [-1,10]
surf_y = [0,0]
elec_x = np.linspace(2,7,5)
elec_y = [0]*5
string1x = [1]*10
string1y = np.linspace(-1,-10,10)
string2x = [8]*10
string2y = np.linspace(-1,-10,10)
poly1x = [3,3,6,6]
poly1y = [-2,-5,-5,-2]
bound1x = [3.5,4.5,5.5]
bound1y = [-6,-7,-8]
buried1x = elec_x
buried1y = np.array(elec_y) - 10
 
geom_input = {'surface': [surf_x,surf_y],
              'electrode':[elec_x,elec_y],
              'borehole1':[string1x,string1y],
              'borehole2':[string2x,string2y],
              'boundary1':[bound1x,bound1y],
              'polygon1':[poly1x,poly1y],
              'buried1':[buried1x,buried1y]}

bh_mesh = mt.tri_mesh(geom_input)

bh_mesh.show()

bh_mesh.write_vtk()
bh_mesh.write_dat()


