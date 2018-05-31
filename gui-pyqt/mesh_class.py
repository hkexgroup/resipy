
"""
Created on Wed May 30 10:19:09 2018, python 3.6.5
Import a vtk file with an unstructured grid (triangular elements) and 
creates a mesh object. 


Functions: 
    tri_cent() - computes the centre point for a 2d triangular element
    vtk_import() - imports a triangular unstructured grid 
Class: 
    mesh_obj
    
@author: jamyd91 (Jimmy Boyd)
"""
#import python standard libraries 
import math as ma
import meshTools as mt
import numpy as np
import time
#import anaconda libraries
from numpy import array
import matplotlib.pyplot as plt
import matplotlib.cm as cmaps
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

#%% create mesh object
class mesh_obj: 
    #create a mesh class
    #put class variables here 
    no_attributes = 1 # it follows we may want to add "attributes to each cell"
    #... we begin assuming each cell has a resistivity assocaited with it but
    #... we may also want associate each cell with a sensitivity for example
    
    def __init__(self,#function constructs our mesh object. 
                 num_nodes,#number of nodes
                 num_elms,#number of elements 
                 node_x,#x coordinates of nodes 
                 node_y,#y coordinates of nodes
                 node_z,#z coordinates of nodes 
                 node_id,#node id number 
                 elm_id,#element id number 
                 node_data,#nodes of element vertices
                 elm_centre,#centre of elements (x,y)
                 elm_area,#area of each element
                 cell_type,#according to vtk format
                 cell_atributes,#the values of the attributes given to each cell 
                 atribute_title,#what is the attribute? we may use conductivity instead of resistivity for example
                 cell_attribute_dump,
                 original_file_path) :
        #assign varaibles to the mesh object 
        self.num_nodes=num_nodes
        self.num_elms=num_elms
        self.node_x=node_x;self.node_y=node_y;self.node_z=node_z
        self.node_id=node_id
        self.elm_id=elm_id
        self.node_data=node_data
        self.elm_centre=elm_centre
        self.elm_area=elm_area
        self.cell_type=cell_type
        self.cell_atributes=cell_atributes 
        self.atribute_title=atribute_title
        #self.cell_attribute_dump=
        self.original_file_path=original_file_path
        
    def file_path(self):#returns the file path from where the mesh was imported
        return(format(self.original_file_path))
       
    def Type2VertsNo(self):#converts vtk cell types into number of vertices each element has 
        if int(self.cell_type[0])==5:#then elements are triangles
            return 3
        elif int(self.cell_type[0])==8 or int(self.cell_type[0])==9:#elements are quads
            return 4
        #add element types as neccessary 
        else:
            print("WARNING: unrecognised cell type")
            return 0
        
    def summary(self):
        #prints summary information about the mesh
        print("\n_______mesh summary_______")
        print("Number of elements: %i"%int(self.num_elms))
        print("Number of nodes: %i"%int(self.num_nodes))
        print("Attribute title: %s"%self.atribute_title)
        print("Number of cell vertices: %i"%self.Type2VertsNo())
        print("Number of cell attributes: %i"%int(self.no_attributes))
        print("original file path: %s"%self.file_path())

    def show(self, ax=None):#displays the mesh using matplotlib
        """
        Show a mesh object using matplotlib. The color map variable should be 
        a string refering to the color map you want (default is "jet").
        As we're using the matplotlib package here any color map avialable within 
        matplotlib package can be used to display the mesh here also. See: 
        https://matplotlib.org/2.0.2/examples/color/colormaps_reference.html
        """ 
        
#        ax.plot(np.random.randn(100,5),'o--')
        
        color_map = 'jet'
        start = time.time()
        if not isinstance(color_map,str):#check the color map variable is a string
            raise NameError('color_map variable is not a string')
            #not currently checking if the passed variable is in the matplotlib library
        
        patches=[]#list wich will hold the polygon instances 
        no_verts=self.Type2VertsNo()#number of vertices each element has 
        for i in range(self.num_elms):
            node_coord=[]#coordinates of the corner of each element 
            for k in range(no_verts):
                node_coord.append((
                        self.node_x[self.node_data[k][i]],
                        self.node_y[self.node_data[k][i]]))                 
            polygon= Polygon(node_coord,True)
            patches.append(polygon) #patch list   
        #build colour map
        X=self.cell_atributes
        colour_array=cmaps.jet(plt.Normalize(min(X),max(X))(X))#maps color onto mesh
        plt.set_cmap(color_map)
        #compile polygons patches into a "patch collection"
        pc=PatchCollection(patches,alpha=0.8,edgecolor='k',facecolor=colour_array)
        pc.set_array(array(X))
        #make figure
#        fig,ax=plt.subplots()#blit polygons to axis
        if ax is None:
            fig, ax = plt.subplots()
        else:
            print('ax is spectified')
#            fig = plt.figure()
        ax.set_aspect('equal')
        ax.add_collection(pc)
        #were dealing with patches and matplotlib isnt smart enough to know what the right limits are 
        ax.set_ylim([min(self.node_y),max(self.node_y)])
        ax.set_xlim([min(self.node_x),max(self.node_x)])
        #update the figure
#        cbar=plt.colorbar(pc,ax=ax)#add colorbar
#        cbar.set_label(self.atribute_title) #set colorbar title      
#        plt.show()
        print('elapsed', time.time()-start)
#        return fig
#        fig.canvas.draw()
        
    def log10(self):#adds a log 10 (resistivity) to the mesh
        #currently this expects that the 
        mesh_obj.no_attributes += 1
        self.log_attribute=[0]*int(self.num_elms)
        for i in range(self.num_elms):
            self.log_attribute[i]=(ma.log10(self.cell_atributes[i]))
            
    def add_attribute(self):
        mesh_obj.no_attributes += 1
        pass #allows us to add an attribute to each element. 

    
    @classmethod # creates a mesh object from a mesh dictionary
    def mesh_dict2obj(cls,mesh_info):
        #check the dictionary is a mesh
        try: 
            if mesh_info['dict_type']!='mesh_info':
                raise NameError("dictionary is not a mesh type")
        except KeyError:
                raise ImportError("dictionary has no dict type variable") 
        #covert into an object 
        obj=cls(mesh_info['num_nodes'],
                     mesh_info['num_elms'], 
                     mesh_info['node_x'],
                     mesh_info['node_y'],
                     mesh_info['node_z'],
                     mesh_info['node_id'],
                     mesh_info['elm_id'],
                     mesh_info['node_data'],
                     mesh_info['elm_centre'],
                     mesh_info['elm_area'],
                     mesh_info['cell_type'],
                     mesh_info['parameters'],
                     mesh_info['parameter_title'],
                     mesh_info['cell_attribute_dump'],
                     mesh_info['original_file_path'])
        return (obj)
    
    @staticmethod
    def help_me():#a basic help me file, needs fleshing out
        available_functions=["show","summary","show_mesh","log10","add_attribute","mesh_dict2obj","Type2VertsNo"]
        print("\n_______________________________________________________")#add some lines, make info look pretty
        print("available functions within the mesh_obj class: \n")
        for i in range(len(available_functions)):
            print("%s"%available_functions[i])
        print("_______________________________________________________")

#%% actual script 
#load in an example mesh 
#mesh_dict=mt.vtk_import(r'C:\Users\blanchy\Downloads\pyScripts\r2gui\f001_res.vtk')#makes a dictionary of a mesh
###
####%% 
#mesh = mesh_obj.mesh_dict2obj(mesh_dict)# this is a mesh_obj class instance 
####print(mesh)#show the mesh object (note objects are not shown in the spyder variable explorer)
###
####%% show some diagonistic information
#mesh.summary()
##mesh.help_me()
#
#fig, ax = plt.subplots()
#mesh.show(ax=ax)
#
#fig.show()
