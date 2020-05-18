# MESH TOOLS 
"""
Created on Wed May 30 10:19:09 2018, python 3.6.5
@author: jamyd91
Module handles mesh generation, display, discretisation and post processing. 
The convention for x y z coordinates is that the z coordinate is the elevation.

Dependencies: 
    numpy (conda lib)
    matplotlib (conda lib)
    gmshWrap(ResIPy resipy module)
    python3 standard libaries
"""
#import standard python packages
import os, platform, warnings, multiprocessing, re, sys
from subprocess import PIPE, Popen
import tempfile
import time, ntpath
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import ListedColormap
import matplotlib.tri as tri
import matplotlib.patches as mpatches
import matplotlib.path as mpath
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import cKDTree
#import R2gui API packages 
import resipy.gmshWrap as gw
from resipy.sliceMesh import sliceMesh # mesh slicing function
import resipy.interpolation as interp

try:#meshCalc needs to be compiled 
    import resipy.cext.meshCalc as mc
except ModuleNotFoundError:
    warnings.warn('meshCalc not installed, meshing options will be limited')

try:#import pyvista if avaiable
    import pyvista as pv
    pyvista_installed = True
except ModuleNotFoundError:
    pyvista_installed = False
    warnings.warn('pyvista not installed, 3D meshing viewing options will be limited')
    
        
#%% cropSurface function
def cropSurface(triang, xsurf, ysurf):
    # check all centroid are below the surface
    trix = np.mean(triang.x[triang.triangles], axis=1)
    triy = np.mean(triang.y[triang.triangles], axis=1)
    
    i2keep = np.ones(len(trix), dtype=bool)
    for i in range(len(xsurf)-1):
        ilateral = (trix > xsurf[i]) & (trix <= xsurf[i+1])
        iabove = (triy > np.min([ysurf[i], ysurf[i+1]]))
        ie = ilateral & iabove
        i2keep[ie] = False
        if np.sum(ie) > 0: # if some triangles are above the min electrode
            slope = (ysurf[i+1]-ysurf[i])/(xsurf[i+1]-xsurf[i])
            offset = ysurf[i] - slope * xsurf[i]
            predy = offset + slope * trix[ie]
            ie2 = triy[ie] < predy # point is above the line joining continuous electrodes
            i2keep[np.where(ie)[0][ie2]] = True
            
    # outside the survey area
    imin = np.argmin(xsurf)
    i2keep[(trix < xsurf[imin]) & (triy > ysurf[imin])] = False
    imax = np.argmax(xsurf)
    i2keep[(trix > xsurf[imax]) & (triy > ysurf[imax])] = False
    
    i2keep1 = i2keep.copy()

    # check all nodes are below the surface
    trix = triang.x
    triy = triang.y
    
    i2keep = np.ones(len(trix), dtype=bool)
    for i in range(len(xsurf)-1):
        ilateral = (trix > xsurf[i]) & (trix <= xsurf[i+1])
        iabove = (triy > np.min([ysurf[i], ysurf[i+1]]))
        ie = ilateral & iabove
        i2keep[ie] = False
        if np.sum(ie) > 0: # if some triangles are above the min electrode
            slope = (ysurf[i+1]-ysurf[i])/(xsurf[i+1]-xsurf[i])
            offset = ysurf[i] - slope * xsurf[i]
            predy = offset + slope * trix[ie]
            ie2 = triy[ie] <= predy # point is above the line joining continuous electrodes
            i2keep[np.where(ie)[0][ie2]] = True
            
    # outside the survey area
    imin = np.argmin(xsurf)
    i2keep[(trix < xsurf[imin]) & (triy > ysurf[imin])] = False
    imax = np.argmax(xsurf)
    i2keep[(trix > xsurf[imax]) & (triy > ysurf[imax])] = False
    
    i2delete = np.where(~i2keep)[0] # get indices
    xi = np.in1d(triang.triangles[:,0], i2delete)
    yi = np.in1d(triang.triangles[:,1], i2delete)
    zi = np.in1d(triang.triangles[:,2], i2delete)
    i2mask = xi | yi | zi
    
    i2keep2 = ~i2mask
    
    return i2keep1 & i2keep2

#%% determine if points are inside cubiod  
def in_box(x,y,z,xmax,xmin,ymax,ymin,zmax,zmin):
    """
    Determine if a point lies inside a bounding volume 
    Parameters
    ----------
    x: array like, float
        x coordinate of query point
    y: array like, float
        y coordinate of query point 
    z: array like, float
        z coordinate of query point
    volume_data: list 
        contains column of poly_data, in the form (polyx, polyy, polyz)
    ray_cast: float, optional
        determines how the far in the x axis a point is ray casted 
    Returns
    ----------
    inside: boolian, numpy array 
        true indexes where point is inside volume
    """
    
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    #first check x coordinates 
    idx_x_in = (x>xmin) & (x<xmax)
    #check y coordinates 
    idx_y_in = (y>ymin) & (y<ymax)
    #check Z coordinates 
    idx_z_in = (z>zmin) & (z<zmax)
    #finally
    idx = (idx_x_in==True) & (idx_y_in==True) & (idx_z_in==True)
    return idx

#%% write descrete points to a vtk file 
def points2vtk (x,y,z,file_name="points.vtk",title='points'):
    """
    Function makes a .vtk file for some xyz coordinates. optional argument
    renames the name of the file (needs file path also) (default is "points.vtk"). 
    title is the name of the vtk file.
            
    Parameters
    ----------
    x : list, tuple, np array
        X coordinates of points.
    y : list, tuple, np array
        Y coordinates of points.
    z : list, tuple, np array
        Z coordinates of points.
    file_name : string, optional
        Path to saved file, defualts to 'points.vtk' in current working directory.
    title : string, optional
        Title of vtk file.
            
    Returns
    -------
    ~.vtk : file
    """
    #error check
    if len(x) != len(y) or len(x) != len(z):
        raise ValueError('mis-match between vector lengths')
    
    fh=open(file_name,'w');#open file handle
    #add header information
    fh.write('# vtk DataFile Version 3.0\n')
    fh.write(title+'\n')
    fh.write('ASCII\n')
    fh.write('DATASET POLYDATA\n')
    #add data
    fh.write('POINTS      %i double\n'%len(x))
    [fh.write('{:<10} {:<10} {:<10}\n'.format(x[i],y[i],z[i])) for i in range(len(x))]
    fh.close()

        
#%% create mesh object
class Mesh:
    """Mesh class.
    
    Parameters
    ----------
    num_nodes : int
        number of nodes
    num_elms : int 
        number of elements 
    node_x : list, 1d numpy array
        x coordinates of nodes 
    node_y : list, 1d numpy array
        coordinates of nodes
    node_z : list, 1d numpy array
        z coordinates of nodes 
    node_id : list
        node id number (ie 1,2,3,4,...)
    elm_id : list
        element id number 
    node_data : list of lists of ints 
        nodes of element vertices in the form [[node1],[node2],[node3],...], each
        node id should be an integer type. 
    elm_centre : list of lists of floats
        centre of elements (x,y)
    elm_area : list 
        area of each element
    cell_type : list of ints
        code referencing cell geometry (e.g. triangle) according to vtk format
    cell_attributes : list of floats
        the values of the attributes given to each cell 
    atribute_title : string 
        what is the attribute? we may use conductivity instead of resistivity for example
    original_file_path : string, optional
        file path to where the mesh file was originally imported
    regions : optional
        element indexes for a material in the mesh (needs further explanation)
        
    Returns
    -------
    Mesh : class
    """

    #%% mesh creation
    def __init__(self,#function constructs our mesh object. 
                 node_x,#x coordinates of nodes 
                 node_y,#y coordinates of nodes
                 node_z,#z coordinates of nodes 
                 node_id,#node id number 
                 elm_id,#element id number 
                 node_data,#nodes of element vertices
                 cell_type,#according to vtk format
                 cell_attributes,#the values of the attributes given to each cell 
                 atribute_title,#what is the attribute? we may use conductivity instead of resistivity for example
                 original_file_path='N/A',
                 order_nodes=True,# order nodes if True, can be computationally expensive 
                 compute_centre=True): # compute cell centres if true, also expensive for big meshes  
        

        #assign variables to the mesh object 
        self.num_nodes = len(node_x)
        self.num_elms = len(node_data[0])
        self.node_x = node_x
        self.node_y = node_y
        self.node_z = node_z
        self.node_id = node_id
        self.elm_id = elm_id
        self.con_matrix = node_data #connection matrix
        self.cell_type = cell_type
        self.cell_attributes = cell_attributes 
        self.atribute_title = atribute_title
        self.original_file_path = original_file_path
        self.surface = None # surface points for cropping the mesh when contouring
        self.e_nodes = None
        self.iremote = None # specify which electrode is remote
        self.cax = None # store mesh.shwo() output for faster mesh.draw()
        self.zone = np.ones(self.num_elms) # by default all in the same zone
        self.elm_centre = None
        self.elm_area = None
        self.attr_cache={'param': np.arange(self.num_elms)+1,
                         'region':cell_attributes,
                         'zone':np.ones(self.num_elms)}# store attributes values per cell
        self.mesh_title = "2D_R2_mesh"
        self.no_attributes = 0
        self.neigh_matrix = None #neighbour matrix, not usually needed unless for 3d tetrahedra problems 
        self.tri_combo = None
        
        #decide on number of dimensions
        if max(node_y) - min(node_y) < 1e-16 or max(node_z) - min(node_z) < 1e-16: # mesh is probably 2D 
            self.ndims=2
        else:
            self.ndims=3
            self.mesh_title = '3D_R3t_mesh'
            
        if self.ndims == 2 and max(node_z) - min(node_z) == 0: # then the mesh is lightly to be a 2D mesh but with elevation in Y axis
            warnings.warn('Y and Z columns in mesh are flipped so that elevation is in the Z axis, to reverse this use mesh.flipYZ()')
            self.flipYZ()
        
        if compute_centre: #compute the centre of mesh cells
            self.cellCentres()
        if order_nodes: # order nodes 
            self.orderNodes()
    
    def copy(self):
        """Return a copy of mesh object. 
        """
        mesh = Mesh(self.node_x,
                    self.node_y,
                    self.node_z,
                    self.node_id,
                    self.elm_id,
                    self.con_matrix,
                    self.cell_type,
                    self.cell_attributes,
                    self.atribute_title,
                    self.original_file_path,
                    order_nodes=False,# should be already ordered
                    compute_centre=False)# should be already computed 
        mesh.elm_centre = self.elm_centre
        mesh.attr_cache = self.attr_cache.copy()
        return mesh 
    
    
    def flipYZ(self):
        """ Make Y column Z column and vice versa, this is useful for when 
        dealing with 2D meshes. Dimensions are modified in place. 
        """
        node_y_cache = np.array(self.node_y)
        node_z_cache = np.array(self.node_z)
        
        self.node_y = list(node_z_cache) # flip the yz columns 
        self.node_z = list(node_y_cache)
        
        self.cellCentres() # recompute cell centres 
                
    
    #%% mesh attribute handling 
    def add_e_nodes(self, e_nodes):
        """Assign node numbers to electrodes. 
        
        Parameters
        ------------
        e_nodes: array like
            array of ints which index the electrode nodes in a mesh
        """
        self.e_nodes = e_nodes
        self.elec_x = np.array(self.node_x)[np.array(e_nodes, dtype=int)]
        if self.ndims==3:
            self.elec_y = np.array(self.node_y)[np.array(e_nodes, dtype=int)]
        else:
            self.elec_y = np.zeros_like(self.elec_x)
        self.elec_z = np.array(self.node_z)[np.array(e_nodes, dtype=int)]
    
    #add some functions to allow adding some extra attributes to mesh 
    def add_sensitivity(self,values):#sensitivity of the mesh
        if len(values)!=self.num_elms:
            raise ValueError("The length of the new attributes array does not match the number of elements in the mesh")
        self.sensitivities = values
        
    def file_path(self):
        """Returns the file path from where the mesh was imported
        """
        return(format(self.original_file_path))
       
    def type2VertsNo(self):
        """Converts vtk cell types into number of vertices each element has
        """
        if int(self.cell_type[0])==5:#then elements are triangles
            return 3
        elif int(self.cell_type[0])==8 or int(self.cell_type[0])==9:#elements are quads
            return 4
        elif int(self.cell_type[0]) == 11: # elements are voxels
            return 8
        elif int(self.cell_type[0]) == 10:# elements are tetrahedra 
            return 4
        elif int(self.cell_type[0]) == 13: # elements are 3d wedges 
            return 6
        #add element types as neccessary 
        else:
            print("WARNING: unrecognised cell type")
            return 0
        
    def summary(self,flag=True):
        """Prints summary information about the mesh
        """
        self.no_attributes = len(self.attr_cache)
        #returns summary information about the mesh, flagto print info, change to return string
        out = "\n_______mesh summary_______\n"
        out += "Number of elements: %i\n"%int(self.num_elms)
        out += "Number of nodes: %i\n"%int(self.num_nodes)
        out += "Number of cell vertices: %i\n"%self.type2VertsNo()
        out += "Number of cell attributes: %i\n"%int(self.no_attributes)
        out += "Dimensions: %i\n"%int(self.ndims)
        out += "original file path: %s\n"%self.file_path()
        if flag:
            print(out)
        else:
            return out

    def __str__(self):
        #returns the summary function if the object is printed using print()
        return self.summary(flag=False) + self.show_avail_attr(flag=False)
            
    def add_attribute(self,values,key):
        """Add a new attribute to mesh. 
        
        Parameters
        ------------
        values: array like
            must have a length which matches the number of elements. Discrete 
            values which map to the elements. 
        key: str
            Name of the attribute, this will be used to reference to the values
            in other mesh functions. 
        """
        if len(values)!=self.num_elms:
            print(len(values),self.num_elms)
            raise ValueError("The length of the new attributes array (%i) does not match the number of elements in the mesh (%i)"%(len(values),self.num_elms))
        self.no_attributes += 1
        try: 
            self.attr_cache[key]=values #allows us to add an attributes to each element.
        except AttributeError:
            self.attr_cache = {}
            self.attr_cache[key]=values #add attribute 
        
    def show_avail_attr(self,flag=True):
        """Show available attributes in mesh.attr_cache. 
        """
        out = '\n______cell attributes_____\n'
        try: 
            for i,key in enumerate(self.attr_cache.keys()):
                out += key + '\n'
        except:
            out += "None\n"
        if flag:
            print(out)
        else:
            return out
    
    def update_attribute(self,new_attributes,new_title='default'):
        """Allows you to reassign the default cell attribute in the mesh object.  
        """
        if len(new_attributes)!=self.num_elms:
            raise ValueError("The length of the new attributes array does not match the number of elements in the mesh")
        self.cell_attributes=new_attributes
        self.atribute_title=str(new_title)
    
    #%% mesh calculations 
    def orderNodes(self,return_count=False):
        """Order mesh nodes in clockwise fashion 
        
        Parameters
        -----------
        return_count:bool, optional
            Function returns the number of reordered elements, default is False. 
        """
        con_mat = self.con_matrix
        con_mata = np.array(self.con_matrix).T.copy() # new object of the connection matrix to reference
        node_x = self.node_x
        node_y = self.node_y
        node_z = self.node_z
        count = 0
        
        if int(self.cell_type[0])==5:#then elements are triangles
            for i in range(self.num_elms):
                n1=(node_x[con_mat[0][i]],node_z[con_mat[0][i]])#define node coordinates
                n2=(node_x[con_mat[1][i]],node_z[con_mat[1][i]])
                n3=(node_x[con_mat[2][i]],node_z[con_mat[2][i]])
                #compute triangle centre
                if interp.ccw(n1,n2,n3)==1:
                    count+=1
                    #points are clockwise and therefore need swapping round
                    con_mat[1][i] = con_mata[i][0]
                    con_mat[0][i] = con_mata[i][1]
                    
        elif int(self.cell_type[0])==8 or int(self.cell_type[0])==9:#elements are quads
            for i in range(self.num_elms):
                n1=(node_x[con_mat[0][i]],node_z[con_mat[0][i]])#define node coordinates
                n2=(node_x[con_mat[1][i]],node_z[con_mat[1][i]])
                n3=(node_x[con_mat[2][i]],node_z[con_mat[2][i]])
                n4=(node_x[con_mat[3][i]],node_z[con_mat[3][i]])
                
                p = np.array((n1,n2,n3,n4)).T
                #compute triangle centre
                order = interp.order_quad(p[0],p[1])
                if list(order) != [0,1,2,3]:
                    count+=1
                    #points are not counter clockwise
                    con_mat[0][i] = con_mata[i][order[0]]
                    con_mat[1][i] = con_mata[i][order[1]]
                    con_mat[2][i] = con_mata[i][order[2]]
                    con_mat[3][i] = con_mata[i][order[3]]
                    
        elif int(self.cell_type[0]) == 11: # elements are voxels
            #print('Node ordering scheme not avialable with this mesh type')
            return
        
        elif int(self.cell_type[0]) == 10:# elements are tetrahedra 
            for i in range(self.num_elms):
                n1=(node_x[con_mat[0][i]],node_y[con_mat[0][i]],node_z[con_mat[0][i]])#define node coordinates
                n2=(node_x[con_mat[1][i]],node_y[con_mat[1][i]],node_z[con_mat[1][i]])
                n3=(node_x[con_mat[2][i]],node_y[con_mat[2][i]],node_z[con_mat[2][i]])
                n4=(node_x[con_mat[3][i]],node_y[con_mat[3][i]],node_z[con_mat[3][i]])
                
                p = np.array((n1,n2,n3,n4)).T
                #compute triangle centre
                if interp.check_tetra(p[0],p[1],p[2]) == 1:
                    count += 1
                    con_mat[1][i] = con_mata[i][0]
                    con_mat[0][i] = con_mata[i][1]
                    
        elif int(self.cell_type[0]) == 13: # elements are 3d wedges 
            for i in range(self.num_elms):
                n1=(node_x[con_mat[0][i]],node_y[con_mat[0][i]],node_z[con_mat[0][i]])#define node coordinates
                n2=(node_x[con_mat[1][i]],node_y[con_mat[1][i]],node_z[con_mat[1][i]])
                n3=(node_x[con_mat[2][i]],node_y[con_mat[2][i]],node_z[con_mat[2][i]])
                # n4=(node_x[con_mat[3][i]],node_y[con_mat[3][i]],node_z[con_mat[3][i]])
                # n5=(node_x[con_mat[4][i]],node_y[con_mat[4][i]],node_z[con_mat[4][i]])
                # n6=(node_x[con_mat[5][i]],node_y[con_mat[5][i]],node_z[con_mat[5][i]])
                #see if top of triangle is counter-clockwise
                if interp.ccw(n1,n2,n3) == 1: #points are clockwise and therefore need swapping round
                    count += 1
                    con_mat[1][i] = con_mata[i][0] # correct the top of the triangle 
                    con_mat[0][i] = con_mata[i][1]
                    con_mat[3][i] = con_mata[i][4] # correct the bottom of the triangle 
                    con_mat[4][i] = con_mata[i][3]

        
        self.con_matrix = con_mat
        
        if return_count:
            return count
        
    def orderElem(self, param=None):
        """Order elements based on the parameter number. Ideally parameters 
        should be concurrent, and fixed elements should go at the base of the 
        connection matrix
        """
        if param is None:
            param = np.array(self.attr_cache['param'])
        else:
            if len(param)!= self.num_elms:
                raise ValueError('The parameter array does not match the number of elements')
                return

        # fixed elements needs to be at the end AND param number continuous
        val = np.sort(np.unique(param)) # unique parameter number
        if np.sum(param == 0) == 0:
            newparam = 1 + np.arange(len(val)) # continuous new param number
        else:
            newparam = np.arange(len(val))
        newval = newparam[np.searchsorted(val, param)] # https://stackoverflow.com/questions/47171356/replace-values-in-numpy-array-based-on-dictionary-and-avoid-overlap-between-new
        param = newval
        order = np.argsort(param)[::-1]
        
        dims = len(self.con_matrix)        
        con_mat_fix = [np.array([0]*self.num_elms)]*dims
        paramFixed = np.array(param)[order]

        for i in range(dims):
            arr = np.array(self.con_matrix[i])
            con_mat_fix[i] = list(arr[order])
            
        for key in self.attr_cache.keys():
            self.attr_cache[key] = list(np.array(self.attr_cache[key])[order])
        
        self.con_matrix = con_mat_fix
        self.cellCentres() # recompute them too
        self.add_attribute(paramFixed,'param')
        
    
    def cellCentres(self):
        """A numpy-based approximation of cell centres for 2D and 3D elements. 
        It's calculated from the mean of cell x y z node coordinates 
        """
        con_mat = np.array(self.con_matrix).T
        vertx = np.array(self.node_x)[con_mat]
        verty = np.array(self.node_y)[con_mat]
        vertz = np.array(self.node_z)[con_mat]

        elm_centre = [[0]*self.num_elms, [0]*self.num_elms, [0]*self.num_elms]
        
        elm_centre[0] = list(np.mean(vertx,axis=1))
        elm_centre[1] = list(np.mean(verty,axis=1))
        elm_centre[2] = list(np.mean(vertz,axis=1))
        self.elm_centre = elm_centre
    
    
    def cellArea(self):
        """Compute the element areas, or in the case of 3D meshes compute the 
        cell volumes. Not yet implimented. 
        """
        con_mat = self.con_matrix
        self.elm_area=[0]*self.num_elms

        if int(self.cell_type[0])==5:#then elements are triangles
            for i in range(self.num_elms):
                n1=(self.node_x[con_mat[0][i]],self.node_z[con_mat[0][i]])#define node coordinates
                n2=(self.node_x[con_mat[1][i]],self.node_z[con_mat[1][i]])
                n3=(self.node_x[con_mat[2][i]],self.node_z[con_mat[2][i]])
                #compute area (for a triangle this is 0.5*base*height)
                base=(((n1[0]-n2[0])**2) + ((n1[1]-n2[1])**2))**0.5
                mid_pt=((n1[0]+n2[0])/2,(n1[1]+n2[1])/2)
                height=(((mid_pt[0]-n3[0])**2) + ((mid_pt[1]-n3[1])**2))**0.5
                self.elm_area[i] = 0.5*base*height
                
        elif int(self.cell_type[0])==8 or int(self.cell_type[0])==9:#elements are quads
            for i in range(self.num_elms):
                n1=(self.node_x[con_mat[0][i]],self.node_z[con_mat[0][i]])#define node coordinates
                n2=(self.node_x[con_mat[1][i]],self.node_z[con_mat[1][i]])
                n3=(self.node_x[con_mat[2][i]],self.node_z[con_mat[2][i]])
                n4=(self.node_x[con_mat[3][i]],self.node_z[con_mat[3][i]])
                p = np.array((n1,n2,n3,n4)).T
                dx = abs(max(p[0]) - min(p[0]))
                dz = abs(max(p[1]) - min(p[1]))
                self.elm_area[i] = dx*dz
                
        elif int(self.cell_type[0]) == 11: # elements are voxels
            for i in range(self.num_elms):
                n1=(self.node_x[con_mat[0][i]],self.node_y[con_mat[0][i]],self.node_z[con_mat[0][i]])#define node coordinates
                n2=(self.node_x[con_mat[1][i]],self.node_y[con_mat[1][i]],self.node_z[con_mat[1][i]])
                n3=(self.node_x[con_mat[2][i]],self.node_y[con_mat[2][i]],self.node_z[con_mat[2][i]])
                n4=(self.node_x[con_mat[3][i]],self.node_y[con_mat[3][i]],self.node_z[con_mat[3][i]])
                n5=(self.node_x[con_mat[4][i]],self.node_y[con_mat[4][i]],self.node_z[con_mat[4][i]])#define node coordinates
                n6=(self.node_x[con_mat[5][i]],self.node_y[con_mat[5][i]],self.node_z[con_mat[5][i]])
                n7=(self.node_x[con_mat[6][i]],self.node_y[con_mat[6][i]],self.node_z[con_mat[6][i]])
                n8=(self.node_x[con_mat[7][i]],self.node_y[con_mat[7][i]],self.node_z[con_mat[7][i]])
                #compute volume (which is a bit of an approximation)
                p = np.array((n1,n2,n3,n4,n5,n6,n7,n8)).T
                dx = abs(max(p[0]) - min(p[0]))
                dy = abs(max(p[1]) - min(p[1]))
                dz = abs(max(p[2]) - min(p[2]))
                self.elm_area[i] = dx*dy*dz

        elif int(self.cell_type[0]) == 10:# elements are tetrahedra 
            print('sorry cant compute element areas for this cell type yet')
            return 
            for i in range(self.num_elms):
                n1=(self.node_x[con_mat[0][i]],self.node_y[con_mat[0][i]],self.node_z[con_mat[0][i]])#define node coordinates
                n2=(self.node_x[con_mat[1][i]],self.node_y[con_mat[1][i]],self.node_z[con_mat[1][i]])
                n3=(self.node_x[con_mat[2][i]],self.node_y[con_mat[2][i]],self.node_z[con_mat[2][i]])
                n4=(self.node_x[con_mat[3][i]],self.node_y[con_mat[3][i]],self.node_z[con_mat[3][i]])

        elif int(self.cell_type[0]) == 13: # elements are 3d wedges 
            for i in range(self.num_elms):
                n1=(self.node_x[con_mat[0][i]],self.node_y[con_mat[0][i]],self.node_z[con_mat[0][i]])#define node coordinates
                n2=(self.node_x[con_mat[1][i]],self.node_y[con_mat[1][i]],self.node_z[con_mat[1][i]])
                n3=(self.node_x[con_mat[2][i]],self.node_y[con_mat[2][i]],self.node_z[con_mat[2][i]])
                n4=(self.node_x[con_mat[3][i]],self.node_y[con_mat[3][i]],self.node_z[con_mat[3][i]])
                n5=(self.node_x[con_mat[3][i]],self.node_y[con_mat[3][i]],self.node_z[con_mat[3][i]])
                n6=(self.node_x[con_mat[3][i]],self.node_y[con_mat[3][i]],self.node_z[con_mat[3][i]])
                #compute wedge volume by computing face area first
                base=(((n1[0]-n2[0])**2) + ((n1[1]-n2[1])**2))**0.5
                mid_pt=((n1[0]+n2[0])/2,(n1[1]+n2[1])/2)
                height=(((mid_pt[0]-n3[0])**2) + ((mid_pt[1]-n3[1])**2))**0.5
                area = 0.5*base*height
                p = np.array((n1,n2,n3,n4,n5,n6)).T
                dz = abs(max(p[2]) - min(p[2]))
                self.elm_area[i] = area * dz

            
    def computeNeigh(self):
        """Compute element neighbour matrix
        """
        if self.type2VertsNo() == 6 or self.type2VertsNo()==8:#not a tetrahedra 3d mesh 
            raise Exception("Sorry neighbour calculation not available yet with this mesh type")
        elif self.type2VertsNo() == 4 and self.ndims==2:#then its a quad mesh 
            raise Exception("Sorry neighbour calculation not available yet with this mesh type")
            
        con_mat = tuple(self.con_matrix)
        if self.ndims == 3:
            self.neigh_matrix, self.tri_combo = mc.neigh3d(con_mat,1)
        elif self.ndims == 2:
            self.neigh_matrix, self.tri_combo = mc.neigh2d(con_mat,1)
        else:
            return 
        
    def refine(self):
        """Refine the mesh into smaller elements 
        """
        print('refining mesh ...', end='')
        if self.ndims==2 and self.type2VertsNo()==3:
            #then its a triangle mesh
            self.splitTri()
        elif self.ndims==2 and self.type2VertsNo()==4:
            # then its a quad mesh (returns a new mesh)
            return self.quad2tri()
        elif self.ndims==3 and self.type2VertsNo()==4:
            #then its a tetra mesh 
            self.splitTetra()
        else:
            print('Sorry not implimented for this mesh type yet')
            return 
        self.cellCentres()
        self.orderNodes()
        print('done')
        
        
    def splitTri(self,param=None):
        """Refine triangles by splitting them into 4 smaller triangles 
        """
        #error checks 
        if self.ndims==3:
            raise ValueError("This kind of mesh splitting isn't available for 3D meshes")
        elif self.type2VertsNo() == 4: #not a triangle mesh 
            raise Exception("Sorry mesh splitting not avialable for this mesh type")
        
        #see if parameter already assigned 
        if param is None:
            if 'param' not in self.attr_cache.keys():
                self.attr_cache['param'] = 1 + np.arange(self.num_elms)
            param = self.attr_cache['param'].copy()
        else:
            if len(param)!= self.num_elms:
                raise ValueError('The parameter array does not match the number of elements')
                return
        
        (new_con_mat, nnode_x, nnode_y, nnode_z,
         nnum_elms, nnum_nodes) = mc.split_tri(
            list(self.con_matrix), 
            list(self.node_x), 
            list(self.node_y), 
            list(self.node_z))
         
        self.con_matrix = new_con_mat
        self.node_x = nnode_x
        self.node_y = nnode_y
        self.node_z = nnode_z
        self.num_elms = nnum_elms
        self.num_nodes = nnum_nodes
        self.cell_attributes = np.repeat(self.cell_attributes,4)
        
        if self.attr_cache is not None:
            for key in self.attr_cache.keys():
                a = np.repeat(self.attr_cache[key],4)
                self.attr_cache[key] = a
        newparam = np.repeat(param,4)

        self.add_attribute(newparam,'param')
             
        
    def splitTetra(self, param=None):
        """Refine tetrahedra by splitting them in six smaller tetrahedra.
        """
        #error checks 
        if self.ndims == 2:
            raise ValueError("This kind of mesh splitting isn't available for 2D meshes")
        if self.type2VertsNo() == 6:#not a tetrahedra 3d mesh 
            raise Exception("This kind of mesh splitting isn't available for meshes of type 'prism'")
            
        #see if parameter already assigned 
        if param is None:
            if 'param' not in self.attr_cache.keys():
                self.attr_cache['param'] = 1 + np.arange(self.num_elms)
            param = self.attr_cache['param'].copy()
        else:
            if len(param)!= self.num_elms:
                raise ValueError('The parameter array does not match the number of elements')
                return
        
        (new_con_mat, nnode_x, nnode_y, nnode_z,
         nnum_elms, nnum_nodes) = mc.split_tetra(
            list(self.con_matrix), 
            list(self.node_x), 
            list(self.node_y), 
            list(self.node_z))
         
        self.con_matrix = new_con_mat
        self.node_x = nnode_x
        self.node_y = nnode_y
        self.node_z = nnode_z
        self.num_elms = nnum_elms
        self.num_nodes = nnum_nodes
        self.cell_attributes = np.repeat(self.cell_attributes,8)
        
        if self.attr_cache is not None:
            for key in self.attr_cache.keys():
                a = np.repeat(self.attr_cache[key],8)
                self.attr_cache[key] = a
        newparam = np.repeat(param,8)
        self.add_attribute(newparam,'param')
    
    def quad2tri(self,param=None):
        """Make a triangle mesh from a quad mesh 

        Returns
        -------
        mesh: class
            Triangle mesh class
        """
        #error checks 
        if self.ndims==3:
            raise ValueError("This kind of mesh splitting isn't available for 3D meshes")
        elif self.type2VertsNo() == 3: #not a triangle mesh 
            raise Exception("Mesh must be quadalateral")
        
        con_mata = np.array(self.con_matrix).T # get connection matrix 
        con1 = con_mata[:,[0,1,2]] # retrieve vertices of columns 
        con2 = con_mata[:,[0,3,2]] # and split to make right angle triangles 
        out_mat = np.zeros((self.num_elms*2,3),dtype=int)       
        out_mat[0::2,:] = con1
        out_mat[1::2,:] = con2 
        
        tri_mesh = self.copy()
        tri_mesh.con_matrix = [list(out_mat[:,0]),list(out_mat[:,1]),list(out_mat[:,2])]
        tri_mesh.num_elms = self.num_elms*2
        tri_mesh.cell_type=[5]
        
        tri_mesh.cell_attributes = np.repeat(self.cell_attributes,2)
        
        tri_mesh.cellCentres()
        
        if self.attr_cache is not None:
            for key in self.attr_cache.keys():
                a = np.repeat(self.attr_cache[key],2)
                tri_mesh.attr_cache[key] = a
                
        if param is None:
            if 'param' not in self.attr_cache.keys():
                self.attr_cache['param'] = 1 + np.arange(self.num_elms)
            param = self.attr_cache['param'].copy()
            
        newparam = np.repeat(param,2)
        tri_mesh.add_attribute(newparam,'param')
        tri_mesh.orderNodes()
        
        return tri_mesh
        
        
    def elemDist(self):
        """Work out the distance of each cell in the mesh to its closest electrode
        """
        elec = np.array([self.elec_x,self.elec_y,self.elec_z]).T
        points = np.array([self.elm_centre[0],
                           self.elm_centre[1],
                           self.elm_centre[1]]).T
        tree = cKDTree(elec) ### setup tree
        dist,idx = tree.query(points) ### >> maps to points to nearest electrode
        self.add_attribute(dist,'cell_distance')
        return dist
    
    def extractSurface(self):
        """ Extract the surface of a triangle or tetrahedral mesh. Ouput of 
        function will depend on mesh type. 
    
        Returns
        -------
        mesh: class
            2d faces of the top of the mesh, if the input mesh is a tetrahedral 
            mesh
        (x,z): tuple
            1D faces of the top of the mesh, if the input is a triangular mesh 

        """
        if self.ndims==2:
            typ = 2
            if self.type2VertsNo() == 4: 
                qmesh = self.quad2tri() # get the triangle mesh instead 
                return qmesh.extractSurface() # run the same function 

        elif self.ndims==3: 
            typ = 3
            if self.type2VertsNo() == 6 or self.type2VertsNo() == 8:#not a tetrahedra 3d mesh 
                raise Exception("Sorry surface extraction isnt available for this type of mesh")
        
        con_mat=self.con_matrix
        if self.neigh_matrix is None: # compute neighbour matrix 
            self.computeNeigh() # this will find the element neighbours and the elements which lie on the outside of the mesh! 
        neigh = np.array(self.neigh_matrix).T
        
        out_elem = np.min(neigh, axis=1) == -1 # elements which have a face on the outside of the mesh 
        neigh_trunc = neigh[out_elem]
        con_trunc = np.array(con_mat).T[out_elem]           
        
        node_x = np.array(self.node_x)
        node_y = np.array(self.node_y)
        node_z = np.array(self.node_z)
        
        if typ == 2: ### Extract 2D faces of top of the mesh ###
            xm = self.elm_centre[0]
            zm = self.elm_centre[2]
            map1 = np.array([1,2]) -1 
            map2 = np.array([2,3]) -1 
            map3 = np.array([3,1]) -1 
            nmap = np.array([map1,map2,map3])
            
            face_list = []
            ocheck = [] # orientation check 
            attr_idx = [] # references attribut index 
            
            x_surf = []
            z_surf = []
            for i in range(len(neigh_trunc)):
                idx = np.argwhere(neigh_trunc[i]==-1)
                for j in range(len(idx)):
                    enodes = con_trunc[i] # element nodes 
                    fnodes = con_trunc[i][nmap[idx[j]]][0] # face nodes
                    face_list.append(fnodes)
                    
                    #now to compute if face boundary condition              
                    xm = np.mean(node_x[fnodes]) #determine approx middle of face
                    zm = np.mean(node_z[fnodes])
                    maxx = np.max(node_x[enodes])
                    maxz = np.max(node_z[enodes])
                    minx = np.min(node_x[enodes])
                    minz = np.min(node_z[enodes])
                    dx = abs(maxx-minx) # work out approx face dimensions 
                    dz = abs(maxz-minz) # work out approx face dimensions 
                    
                    if dx < 1e-16: # element on side of mesh 
                        ocheck.append(0) 
                    else: 
                        x = np.append(node_x[fnodes],xm)
                        z = np.append(node_z[fnodes],zm+dx+dz)
                        o = interp.ccw([x[0],z[0]], 
                                       [x[1],z[1]], 
                                       [x[2],z[2]])
                        nacheckidx = (node_x >= minx) & (node_x <= maxx) # node above check index 
                        nacheck = node_z[nacheckidx] > maxz
                        # return nacheck,node_z[nacheckidx] > zm, nacheckidx
                        if o == 1 and not any(nacheck==True):
                            x_surf.append(x[0])
                            x_surf.append(x[1])
                            z_surf.append(z[0])
                            z_surf.append(z[1])     
                        ocheck.append(o)
                    attr_idx.append(i)

            uni_x, idx, counts = np.unique(x_surf,return_counts=True,return_index=True)
            uni_z = np.array(z_surf)[idx]
            return (uni_x,uni_z)
            
        else:  ### Extract 3D faces of top of the mesh ###
            ## map used in neigh calculation ##         
            map1 = np.array([2,3,4]) -1 
            map2 = np.array([1,4,3]) -1 
            map3 = np.array([1,2,4]) -1 
            map4 = np.array([1,2,3]) -1 
            nmap = np.array([map1,map2,map3,map4])
            
            face_list = []
            ocheck = [] # orientation check 
            attr_idx = [] # references attribut index 
            for i in range(len(neigh_trunc)):
                idx = np.argwhere(neigh_trunc[i]==-1)
                for j in range(len(idx)):
                    enodes = con_trunc[i] # element nodes 
                    fnodes = con_trunc[i][nmap[idx[j]]][0] # face nodes
                    mnode = [enodes[k] not in fnodes for k in range(4)] # find node missing from the face 
                    x = np.append(node_x[fnodes],node_x[enodes[mnode]])
                    y = np.append(node_y[fnodes],node_y[enodes[mnode]])
                    z = np.append(node_z[fnodes],node_z[enodes[mnode]])
                    if interp.check_tetra(x,y,z)==2:#points are counter clockwise
                        fnodes_sorted = fnodes
                    else: # reorganise so they are counter clockwise
                        fnodes_sorted = np.array((fnodes[2],fnodes[1],fnodes[0]))
                    face_list.append(fnodes_sorted)
                    
                    #now to compute if face boundary condition              
                    xm = np.mean(node_x[fnodes]) #determine approx middle of face
                    ym = np.mean(node_y[fnodes])
                    zm = np.mean(node_z[fnodes])
                    dx = abs(max(node_x[fnodes]) - min(node_x[fnodes])) # work out approx face dimensions 
                    dy = abs(max(node_y[fnodes]) - min(node_y[fnodes]))
                    dz = abs(max(node_z[fnodes]) - min(node_z[fnodes]))
                    
                    if dy < 1e-16 or dx < 1e-16 :
                        ocheck.append(0) 
                    else: 
                        x = np.append(node_x[fnodes_sorted],xm)
                        y = np.append(node_y[fnodes_sorted],ym)
                        z = np.append(node_z[fnodes_sorted],zm+dx+dy+dz)
                        o = interp.check_tetra(x,y,z)

                        if o == 1:# check if nodes exit exist above the element
                            path = mpath.Path(np.array([node_x[fnodes],node_y[fnodes]]).T)
                            inside = path.contains_points(np.array([node_x,node_y]).T) 
                            nacheck = node_z[inside] > max(node_z[fnodes]) 
                            if any(nacheck==True):
                                o = 2
                        
                        ocheck.append(o)
                    #if ocheck == 0 then the face is on the side of the mesh 
                    #if ocheck == 1 then the face is on top of the mesh 
                    #if ocheck == 2 then the face is on the bottom of the mesh (or inside the mesh) 
                    attr_idx.append(i)
            
            #assign node boundary markers 
            ochecka = np.array(ocheck,dtype=int)
            ikeep = ochecka == 1
            face_matrix = np.array(face_list)[ikeep,:]
            attr_idx = np.array(attr_idx,dtype=int)[ikeep]
            
            face_cmat = [list(face_matrix[:,0]),
                         list(face_matrix[:,1]),
                         list(face_matrix[:,2])]
            no_elms = len(face_matrix) # number of elements 

            nmesh = Mesh(node_x, # make new mesh 
                         node_y, 
                         node_z, 
                         node_id = np.arange(len(node_x)), 
                         elm_id = np.arange(no_elms), 
                         node_data = face_cmat, 
                         cell_type = [5], 
                         cell_attributes = [0]*no_elms, 
                         atribute_title='face')
            
            for key in self.attr_cache.keys(): # add attributes 
                X = np.array(self.attr_cache[key])[out_elem]
                nmesh.add_attribute(X[attr_idx], key)                
            nmesh.__rmexcessNodes() # remove excess nodes which are not used 
            
            return nmesh
        
    #%% Truncating the mesh 
    
    def __rmexcessNodes(self):
        """ Remove any nodes are not inside the connection matrix
        """
        con_mat = np.array(self.con_matrix).T.copy()
        shp = con_mat.shape # connection matrix shape 
        con_matf = con_mat.flatten() # flattened connection matrix 
        sort_idx = np.argsort(con_matf) # sort the connection matrix 
        map_idx = np.argsort(sort_idx)# indexes to map the indexes back to thier original positions 
        sortf = con_matf[sort_idx] # sorted connection matrix 
        uni_nodes = np.unique(con_matf) # unique mesh nodes 
        
        #create new node indexes
        new_nodes = np.array([0]*len(sortf))
        count = 0
        new_nodes[0] = count
        for i in range(1,len(sortf)):
            if sortf[i] != sortf[i-1]:
                count += 1
            new_nodes[i] = count
        
        #remap indexes 
        self.node_x = list(np.array(self.node_x)[uni_nodes])
        self.node_y = list(np.array(self.node_y)[uni_nodes])
        self.node_z = list(np.array(self.node_z)[uni_nodes])
        new_conf = new_nodes[map_idx]
        new_con = new_conf.reshape(shp)
        for i in range(shp[1]):
            self.con_matrix[i] = list(new_con[:,i])
            
        self.num_nodes = len(uni_nodes)
        
    def crop(self, polyline):
        """Crop the mesh given a polyline in 2D.
        
        Parameters
        ----------
        polyline : array of float
            Array of size Nx2 with the XZ coordinates forming the polyline. Note
            that the first and last coordinates should be the
            same to close the polyline.
        """
        # get points inside the polygon
        path = mpath.Path(polyline)
        centroids = np.c_[self.elm_centre[0], self.elm_centre[2]]
        i2keep = path.contains_points(centroids) # this is the fastest way to check if point are in polygon 
        
        # filter element-based attribute
        ecx = np.array(self.elm_centre[0])[i2keep].tolist()
        ecy = np.array(self.elm_centre[1])[i2keep].tolist()
        ecz = np.array(self.elm_centre[2])[i2keep].tolist()
        self.elm_centre = (ecx, ecy, ecz)
        self.con_matrix = tuple(np.array(self.con_matrix)[:,i2keep].tolist())
        self.elm_area = np.array(self.elm_area)[i2keep].tolist()
        self.elm_id = np.array(self.elm_id)[i2keep].tolist()
        self.cell_attributes = np.array(self.cell_attributes)[i2keep].tolist()
        for key in self.attr_cache.keys():
            self.attr_cache[key] = np.array(self.attr_cache[key])[i2keep]
        self.num_elms = np.sum(i2keep)
        
        self.__rmexcessNodes            
        
        # return a truncated mesh 
    def truncateMesh(self,xlim=None,ylim=None,zlim=None):
        """Crop the mesh to a box of given limits, like how R3t behaves 
        when outputting inverted results. 
        
        Parameters
        ------------
        xlim : tuple, optional
            Axis x limits as `(xmin, xmax)`.
        ylim : tuple, optional
            Axis y limits as `(ymin, ymax)`. 
        zlim : tuple, optional
            Axis z limits as `(ymin, ymax)`. 
            
        Returns
        ------------
        mesh: Class
            New instance of mesh class which is truncated 
        """
        if xlim is None:
            xlim=[min(self.elec_x), max(self.elec_x)]
        if ylim is None:
            ylim=[min(self.elec_y), max(self.elec_y)]
        if zlim is None:
            zlim=[min(self.elec_z), max(self.elec_z)]
            
        elm_x = np.array(self.elm_centre[0])
        elm_y = np.array(self.elm_centre[1])
        elm_z = np.array(self.elm_centre[2])
        in_elem = in_box(elm_x,elm_y,elm_z,xlim[1],xlim[0],ylim[1],ylim[0],zlim[1],zlim[0])#find inside of limits 
        temp_con_mat = np.array(self.con_matrix,dtype='int64')#temporary connection matrix which is just the elements inside the box
        #con_mat=list(temp_con_mat[:,in_elem]) # truncate connection matrix
        
        new_attr_cache = self.attr_cache.copy()
        new_attr = np.array(self.cell_attributes)
        
        elm_id = np.array(self.elm_id)

        #truncate the attribute table down to the inside elements 
        for key in self.attr_cache.keys():
            X = np.array(self.attr_cache[key])
            new_attr_cache[key] = X[in_elem]
        
        nmesh = self.copy() # make a new mesh object with fewer elements 
        
        nmesh.attr_cache = new_attr_cache
        nmesh.num_elms = len(elm_id[in_elem])
        nmesh.cell_attributes = new_attr[in_elem]
        nmesh.elm_id = elm_id[in_elem]
        new_elm_centre = [[],[],[]]
        new_elm_centre[0] = elm_x[in_elem]
        new_elm_centre[1] = elm_y[in_elem]
        new_elm_centre[2] = elm_z[in_elem]
        nmesh.elm_centre = new_elm_centre
        
        new_con_mat = temp_con_mat[:,in_elem].copy()            
        nmesh.con_matrix=list(new_con_mat)
        
        nmesh.__rmexcessNodes() # remove the excess nodes 
            
        return nmesh # return truncated mesh 
    
    def threshold(self,attr=None,vmin=None,vmax=None):
        """Threshold the mesh to certian attribute values. 
        
        Parameters
        ------------
        attr: string
            Name of attribute to threshold by 
        vmin: float
            minimum value of attribute
        vmax: float
            maximum value of attribute
        
        Returns
        ------------
        mesh: Class
            New instance of mesh class which is thresholded 
        """
        if self.attr_cache is None:
            raise EnvironmentError("No 'attr_cache' varaible exists for the mesh class!")
        elif attr not in self.attr_cache.keys():
            raise ValueError("Specified attribute has not been defined")
            
        X = np.array(self.attr_cache[attr])
        
        if vmin is None:
            vmin = np.min(X)
        if vmax is None:
            vmax = np.max(X)
            
        in_elem = (X >= vmin) & (X <= vmax) # index of elements to keep 
        
        temp_con_mat = np.array(self.con_matrix,dtype='int64')#temporary connection matrix which is just the elements inside the box
        
        new_attr_cache = self.attr_cache.copy()
        new_attr = np.array(self.cell_attributes)
        
        elm_id = np.array(self.elm_id)

        #truncate the attribute table down to the inside elements 
        for key in self.attr_cache.keys():
            X = np.array(self.attr_cache[key])
            new_attr_cache[key] = X[in_elem]
        
        nmesh = self.copy() # make a new mesh object with fewer elements 
        
        nmesh.attr_cache = new_attr_cache
        nmesh.num_elms = len(elm_id[in_elem])
        nmesh.cell_attributes = new_attr[in_elem]
        nmesh.elm_id = elm_id[in_elem]
        new_elm_centre = [[],[],[]]
        new_elm_centre[0] = list(np.array(self.elm_centre[0])[in_elem])
        new_elm_centre[1] = list(np.array(self.elm_centre[1])[in_elem])
        new_elm_centre[2] = list(np.array(self.elm_centre[2])[in_elem])
        nmesh.elm_centre = new_elm_centre
        
        new_con_mat = temp_con_mat[:,in_elem].copy()            
        nmesh.con_matrix=list(new_con_mat)
        
        nmesh.__rmexcessNodes() # remove the excess nodes 
            
        return nmesh # return truncated mesh    

#%% lookup information      
    def quadMeshNp(self, topo=None):
        """ Convert mesh nodes into x column indexes in the case of quad meshes. 
        Does not currently support changes in electrode elevation! 
        
        Returns
        ----------
        colx: list
            X column indexes for quad mesh 
        """
        if int(self.cell_type[0])!=8 and int(self.cell_type[0])!=9:#elements are not quads
            raise TypeError('Mesh is not composed of 2D quads')
        
        unix = np.unique(self.node_x) # unique x values in the x node coordinates 
        if topo is None: # find the y column is a little challanging without knowing the original topography 
            uniz = np.unique(self.node_z) # if you dont know any better then just use the unique z values 
        else:
            uniz = topo # ideally use mesh topography 
        e_nodes = self.e_nodes
        colx = [0]*len(e_nodes) # column indexes for x coordinates 
        colz = [0]*len(e_nodes) # column indexes for z coordinates 
        for i in range(len(e_nodes)):
            x = self.node_x[e_nodes[i]] # get node x coordinate 
            z = self.node_z[e_nodes[i]]
            colx[i] = int(np.argwhere(x==unix) + 1) # find its index 
            colz[i] = int(np.argwhere(z==uniz) + 1)
            
        return colx#,colz # return columns to go in parameters 
        
    def meshLookUp(self,look_up_mesh):
        """Look up values from another mesh using nearest neighbour look up, 
        assign attributes to the current mesh class. 
        
        Parameters
        -------------
        look_up_mesh: class
            Another mesh class. 
        
        Notes
        -------------
        This can fail for  large meshes due to the size of matrices involved. 
        """
        #assign coordinate arrays 
        attr_cache = {} # pd.DataFrame()
        x_old = look_up_mesh.elm_centre[0]
        x_new = self.elm_centre[0]
        y_old = look_up_mesh.elm_centre[1]
        y_new = self.elm_centre[1]
        z_old = look_up_mesh.elm_centre[2]
        z_new = self.elm_centre[2]
        i_old = np.array(look_up_mesh.cell_attributes)
        #do look up 
        if self.ndims==3:
            i_new, idxes = interp.nearest3d(x_new,y_new,z_new,
                                            x_old,y_old,z_old,i_old,
                                            return_idx=True)
        elif self.ndims==2:
            i_new, idxes = interp.nearest(x_new,z_new,
                                            x_old,z_old,i_old,
                                            return_idx=True) 
            
        look_up_cache = look_up_mesh.attr_cache.copy()
        look_up_keys = look_up_cache.keys()
        for k in look_up_keys:
            look_up_array = np.array(look_up_cache[k])
            attr_cache[k] = look_up_array[idxes]
        self.attr_cache = attr_cache.copy()
            
    def trans_mesh(self,x,y,z):
        """Translate mesh by x y z coordinate
        """
        self.node_x = np.array(self.node_x)+x
        self.node_y = np.array(self.node_y)+y
        self.node_z = np.array(self.node_z)+z
        
    #%% mesh display 
    def _clipContour(self, ax, cont, maxDepth=None):
        """Clip contours using mesh bound and surface if available.
        
        Parameters
        ----------
        ax : matplotlib.Axes
            Axis.
        cont : matplotlib.collections
            Collection of contours.
        maxDepth : float (m), optional
                Depth of the fine/coarse region in mesh. 
                If not None: Contour plots will crop out below maxDepth [m].
        """
        # mask outer region
        xmin = np.min(self.node_x)
        xmax = np.max(self.node_x)
        zmin = np.min(self.node_z)
        zmax = np.max(self.node_z)
        if self.surface is not None:
            xsurf, zsurf = self.surface[:,0], self.surface[:,1]
            if maxDepth is not None:
                xfmd, zfmd = self.surface[:,0][::-1], self.surface[:,1][::-1] - maxDepth
                verts = np.c_[np.r_[xmin, xmin, xsurf, xmax, xmax, xfmd, xmin],
                              np.r_[zmin, zmax, zsurf, zmax, zmin, zfmd, zmin]]
            else:
                verts = np.c_[np.r_[xmin, xmin, xsurf, xmax, xmax, xmin],
                              np.r_[zmin, zmax, zsurf, zmax, zmin, zmin]]     
        else:
            verts = np.c_[np.r_[xmin, xmin, xmax, xmax, xmin],
                          np.r_[zmin, zmax, zmax, zmin, zmin]]                
        # cliping using a patch (https://stackoverflow.com/questions/25688573/matplotlib-set-clip-path-requires-patch-to-be-plotted)
        path = mpath.Path(verts)
        patch = mpatches.PathPatch(path, facecolor='none', edgecolor='none')
        ax.add_patch(patch) # need to add so it knows the transform
        for col in cont.collections:
            col.set_clip_path(patch)        

    def show(self,color_map = 'Spectral',#displays the mesh using matplotlib
             color_bar = True,
             xlim = None,
             zlim = None,
             ax = None,
             electrodes = True,
             sens = False,
             edge_color = 'k',
             contour = False,
             vmin = None,
             vmax = None,
             attr = None,
             clabel = None,
             hor_cbar = False,
             sensPrc = None,
             maxDepth = None,
             aspect = 'equal',
             **kwargs):
        """ Displays a 2d mesh and attribute.
        
        Parameters
        ----------
        color_map : string, optional
            color map reference 
        color_bar : Boolean, optional 
            `True` to plot colorbar 
        xlim : tuple, optional
            Axis x limits as `(xmin, xmax)`.
        zlim : tuple, optional
            Axis z limits as `(zmin, zmax)`. 
        ax : matplotlib axis handle, optional
            Axis handle if preexisting (error will thrown up if not) figure is to be cast to.
        electrodes : boolean, optional
            Enter true to add electrodes to plot.
        sens : boolean, optional
            Enter true to plot sensitivities. 
        edge_color : string, optional
            Color of the cell edges, set to `None` if you dont want an edge.
        contour : boolean, optional
            If `True`, plot filled with contours instead of the mesh.
        vmin : float, optional
            Minimum limit for the color bar scale.
        vmax : float, optional
            Maximum limit for the color bar scale.
        attr : string, optional
            Which attribute in the mesh to plot, references a dictionary of attributes. attr is passed 
            as the key for this dictionary.
        clabel : string, optional
            Label of the colorbar. Default is the value of `attr` argument.
        hor_cbar : boolean, optional
            'True' to make a horizontal color bar at the bottom of the plot, default
            is vertical color bar to the right of the plot. 
        sensPrc : float, optional
            Normalised (between 0 and 1) sensitivity value threshold. Default
            is None meaning the sensitivity is just overlay. Need `sens=True` 
            to be used.
        maxDepth : float 
            Maximum absolute depth to be shown on the plotted figure. 
        aspect : string, optional
            defines the aspect ratio of the plot.
            'equal' locks the aspect ratio.
            'auto', aspect ratio is define by plotting area.
        
        Returns
        -------
        figure : matplotlib figure 
            Figure handle for the plotted mesh object.
        
        Notes
        -----
        Show a mesh object using matplotlib. The color map variable should be 
        a string refering to the color map you want (default is "jet").
        As we're using the matplotlib package here any color map avialable within 
        matplotlib package can be used to display the mesh here also. See: 
        https://matplotlib.org/2.0.2/examples/color/colormaps_reference.html
        """
        #check color map argument is a string 
        if not isinstance(color_map,str):#check the color map variable is a string
            raise NameError('color_map variable is not a string')
            #not currently checking if the passed variable is in the matplotlib library
        
        if self.iremote is None:
            try:
                iremote = np.zeros(len(self.elec_x), dtype=bool)
                self.iremote = iremote
            except Exception as e:
                print('No electrode found: ', e)
                pass
        else:
            iremote = self.iremote

        if self.ndims == 3:
            self.show3D(color_map = color_map,#displays the mesh using matplotlib
             color_bar = color_bar, # pass arguments to 3D function
             xlim = xlim,
             zlim = zlim,
             ax = ax,
             electrodes = electrodes,
             sens = sens,
             edge_color = edge_color,
             vmin=vmin,
             vmax=vmax,
             attr=attr,
             **kwargs) # show 3D mesh instead 
            return # exit 2D mesh show function 
            
        
        # decide which attribute to plot, we may decide to have other attritbutes! 
        if attr not in self.attr_cache.keys():
            attr == None
        if attr == None:
            attr = list(self.attr_cache.keys())[3]
            X = np.array(self.attr_cache[attr])
            color_bar_title = attr
        else:
            X = np.array(self.attr_cache[attr])
            color_bar_title = attr
            
        if clabel is not None:
            color_bar_title = clabel

        iplot = False
        if ax is None:
            iplot = True
            fig,ax=plt.subplots()
            self.fig = fig
            self.ax = ax
        else:
            self.fig = ax.figure
            self.ax = ax
        #if no dimensions are given then set the plot limits to edge of mesh
        
        try: 
            elec_x = self.elec_x[~iremote]
            if xlim==None:
                xlim=[min(elec_x),max(elec_x)]
            if zlim==None:
                doiEstimate = 2/3*np.abs(elec_x[0]-elec_x[-1])
                # longest dipole calculation available in R2 class
                zlim=[min(self.elec_z)-doiEstimate,max(self.elec_z)]
        except AttributeError:
            if xlim==None:
                xlim=[min(self.node_x),max(self.node_x)]
            if zlim==None:
                zlim=[min(self.node_z),max(self.node_z)]

        if np.diff(xlim) == 0: # protection against thin axis margins 
            xlim=[xlim[0]-2,xlim[1]+2]
        if np.diff(zlim) == 0:
            zlim=[zlim[0]-2,zlim[1]+2]
                
        ##plot mesh! ##
        #compile mesh coordinates into polygon coordinates  
        nodes = np.c_[self.node_x, self.node_z]
        connection = np.array(self.con_matrix).T # connection matrix 
        if maxDepth is not None: 
            depths = np.array(self.computeElmDepth())
            ikeep = depths < maxDepth
            #truncate connection matrix and plotted array
            connection = connection[ikeep,:]
            X = X[ikeep]
            
        #compile polygons patches into a "patch collection"
        ###X=np.array(self.cell_attributes) # maps resistivity values on the color map### <-- disabled 
        coordinates = nodes[connection]
        if vmin is None:
            vmin = np.min(X)
        if vmax is None:
            vmax = np.max(X)
        
        if edge_color == None or edge_color=='none' or edge_color=='None':
            edge_color='face'#set the edge colours to the colours of the polygon patches

        if contour is False:
            if attr is 'region': # so the default material
                cm = plt.get_cmap(color_map, len(np.unique(X))) # this makes a discrete colormap
            else:
                cm = color_map
            coll = PolyCollection(coordinates, array=X, cmap=cm, edgecolors=edge_color,linewidth=0.5)
            coll.set_clim(vmin=vmin, vmax=vmax)
            ax.add_collection(coll)#blit polygons to axis
#            triang = tri.Triangulation(nodes[:,0], nodes[:,1], connection)
#            coll = ax.tripcolor(triang, X, cmap=color_map, edgecolors=edge_color, linewidth=0.5)
            self.cax = coll

        else:#use contour algorithm (only for 2D and y is considered depth here)
            if maxDepth is not None:
                xc = np.array(self.elm_centre[0])[ikeep]
                yc = np.array(self.elm_centre[2])[ikeep]
            else:
                xc = np.array(self.elm_centre[0])
                yc = np.array(self.elm_centre[2])
            zc = np.array(X)
            
            # check for 0 in sigma log
            if attr == 'Sigma_imag(log10)':
                ie = zc != 0
                zc = zc[ie]
                xc = xc[ie]
                yc = yc[ie]
            x = np.array(self.node_x)
            y = np.array(self.node_z)
            
            # set scale arrangement
            if vmin is None:
                vmin = np.nanmin(zc)
            if vmax is None:
                vmax = np.nanmax(zc)
            if vmax > vmin:
                levels = np.linspace(vmin, vmax, 13) # to have 2 contours between two cbar values!
            else:
                levels = None
            
            if self.cell_type[0] == -1:
                print('really ?!')
                pass
#            if self.cell_type[0] == 9: # quadrilateral mesh (exact topo)   
#                # interpolate the cell-centered value to the node to be able
#                # to use the quadrilateral mesh already in the grid
#                z = interp.linear(x, y, xc, yc, zc)
#                
#                def rebuildRegularGrid(x2, y2, z2):
#                    x2unique = np.unique(x2)
#                    xs = []
#                    for xuni in x2unique: # for loop otherwise dataframe groupby
#                        xs.append(np.where(x2 == xuni)[0])
#                    minLength = np.min([len(a) for a in xs])
#                    xs2 = []
#                    for a in xs:
#                        isort = np.argsort(y2[a])
#                        xs2.append(a[isort][-minLength:])
#                    xs2 = np.vstack(xs2)
#                    X = x2[xs2]
#                    Y = y2[xs2]
#                    Z = z2[xs2]
#                    return X, Y, Z
#            
#                Xi, Yi, Zi = rebuildRegularGrid(x, y, z)
#                self.cax = ax.contourf(Xi, Yi, Zi, levels=levels, cmap=color_map, extend = 'both')
#            
            
#            elif self.cell_type[0] == 5: # triangular mesh (exact topo)
#                # interpolate the cell-centered value to the node to be able
#                # to use the triangular mesh already in the grid
#                z = interp.nearest(x, y, xc, yc, zc)
#                triang = tri.Triangulation(x, y, connection)
#                
#                self.cax = ax.tricontourf(triang, z, levels=levels, extend='both', cmap=color_map)
            
            else: # fallback mode with tricontourf and cropSurface() (topo based on centroids) 
#                if self.surface is not None:
#                    xf, yf = self.surface[:,0], self.surface[:,1]
#                    zf = interp.nearest(xf, yf, xc, yc, zc) # interpolate before overiding xc and yc
#                    xc = np.r_[xc, xf]
#                    yc = np.r_[yc, yf]
#                    zc = np.r_[zc, zf]
#                triang = tri.Triangulation(xc, yc) # build grid based on centroids
#                
#                if self.surface is not None:
#                    try:
#                        triang.set_mask(~cropSurface(triang, self.surface[:,0], self.surface[:,1]))
#                    except Exception as e:
#                        print('Error in Mesh.show() for contouring: ', e)
#                
                triang = tri.Triangulation(xc, yc) # build grid based on centroids
                self.cax = ax.tricontourf(triang, zc, levels=levels, extend='both', cmap=color_map)
                self._clipContour(ax, self.cax, maxDepth=maxDepth)
            
        ax.autoscale()
        #were dealing with patches and matplotlib isnt smart enough to know what the right limits are, hence set axis limits 
        ax.set_ylim(zlim)
        ax.set_xlim(xlim)
        ax.set_xlabel('Distance [m]')
        ax.set_ylabel('Elevation [m]')

        if color_bar:#add the color bar 
            cbar_horizontal = 'vertical'
            if hor_cbar: # change orientation if true 
                cbar_horizontal = 'horizontal'
            self.cbar = plt.colorbar(self.cax, ax=ax, format='%.1f',orientation=cbar_horizontal, fraction=0.046, pad=0.04)
            if attr == 'region': # default to material
                val = np.sort(np.unique(X)).astype(int)
                if len(val) > 1:
                    interval = (val[-1]-val[0])/len(val)
                    self.cbar.set_ticks(np.arange(val[0]+interval/2, val[-1], interval))
                else:
                    self.cbar.set_ticks([1])
                self.cbar.set_ticklabels(val)
            self.cbar.set_label(color_bar_title) #set colorbar title

        ax.set_aspect(aspect)#set aspect ratio equal (stops a funny looking mesh)

        #biuld alpha channel if we have sensitivities 
        if sens:
            if 'Sensitivity(log10)' not in self.attr_cache.keys():
                print('ERROR: No sensitivity attribute found')
            else:
                try:
                    weights = np.array(self.attr_cache['Sensitivity(log10)']) #values assigned to alpha channels 
                    if maxDepth is not None:
                        weights = weights[ikeep]
                    if sensPrc is None:
                        thresh = np.log10(0.001*(10**np.nanmax(weights)))
        #                    thresh = np.percentile(weights, 50, interpolation='nearest')
                        x = np.sort(weights)
                        i = np.where(x > thresh)[0][0]
                        x = np.argsort(weights)
                        alphas = np.zeros(self.num_elms)
                        alphas[:i] = np.linspace(1, 0, len(alphas[:i]))
                        raw_alpha = np.ones((self.num_elms,4),dtype=float) #raw alpha values 
                        raw_alpha[:, -1] = alphas
                        alpha_map = ListedColormap(raw_alpha) # make a alpha color map which can be called by matplotlib
                        #make alpha collection
                        alpha_coll = PolyCollection(coordinates, array=weights, cmap=alpha_map, edgecolors='none', linewidths=0)#'face')
                        #*** the above line can cuase issues "attribute error" no np.array has not attribute get_transform, 
                        #*** i still cant figure out why this is because its the same code used to plot the resistivities 
                        ax.add_collection(alpha_coll)
                    else:
                        #values assigned to alpha channels 
                        a = np.log10(0.000001*(10**np.nanmax(weights)))
                        b = np.log10(0.1*(10**np.nanmax(weights)))
                        ab = np.linspace(a, b, 100)
                        thresh = ab[int(sensPrc*99)]
        #                    thresh = np.percentile(weights, sensPrc*100, interpolation='nearest')
                        x = np.sort(weights)
                        i = np.where(x > thresh)[0][0]
                        x = np.argsort(weights)
                        alphas = np.zeros(self.num_elms)
                        alphas[:i] = np.linspace(1, 0.2, len(alphas[:i]))
                        raw_alpha = np.ones((self.num_elms,4),dtype=float) #raw alpha values 
                        raw_alpha[:, -1] = alphas
                        alpha_map = ListedColormap(raw_alpha) # make a alpha color map which can be called by matplotlib
                        #make alpha collection
                        alpha_coll = PolyCollection(coordinates, array=weights, cmap=alpha_map, edgecolors='none', linewidths=0)#'face')
                        #*** the above line can cuase issues "attribute error" no np.array has not attribute get_transform, 
                        #*** i still cant figure out why this is because its the same code used to plot the resistivities 
                        ax.add_collection(alpha_coll)
                    
                    # legacy implementation of the mesh cropping below (don't delete it)
#                    x = np.array(self.attr_cache['Sensitivity(log10)'])
#                    alphas = np.zeros(len(x))
#                    xmin, xmax = np.nanmin(x), np.nanmax(x)
#                    xnorm = (x-xmin)/(xmax-xmin) # should all be between 0 and 1
#                    i2mask = xnorm < sensPrc
#                    xnorm[i2mask] = 0
#                    xnorm[~i2mask] = 1
#                    raw_alpha = np.ones((2, 4), dtype=float)
#                    raw_alpha[-1,-1] = 0 # larger sensitivity values are transparent
#                    alpha_map = ListedColormap(raw_alpha)
#                    if contour is False:
#                        alpha_coll = PolyCollection(coordinates[i2mask], array=xnorm[i2mask], cmap=alpha_map, edgecolors='face', linewidths=1)
#                        ax.add_collection(alpha_coll)
#                    else: # if contour is True, cropSurface() should have been defined
#                        xc = np.array(self.elm_centre[0])
#                        yc = np.array(self.elm_centre[2])
#                        x = np.array(self.node_x)
#                        y = np.array(self.node_z)
#                        zc = xnorm # normalized sensitivity here
#                        
#                        # doesn't work as the QHull alright discard the elements forming the topo :/
#    #                    from scipy.interpolate import LinearNDInterpolator
#    #                    lin = LinearNDInterpolator(np.c_[xc, yc], zc)
#    #                    z = lin(np.c_[x, y])
#    #                    print(np.sum(np.isnan(z)))
#    #                    inan = ~np.isnan(z)
#    #                    x, y, z = x[inan], y[inan], z[inan]
#    #                    triang = tri.Triangulation(x, y)
#                        
#                        # adding surface points to form surface triangles
#                        if self.surface is not None:
#                            xf, yf = self.surface[:,0], self.surface[:,1]
#                            zf = interp.nearest(xf, yf, xc, yc, zc) # interpolate before overiding xc and yc
#                            xc = np.r_[xc, xf]
#                            yc = np.r_[yc, yf]
#                            zc = np.r_[zc, zf]
#                        triang = tri.Triangulation(xc, yc) # build grid based on centroids and surface points
#                        z = zc
#                    
#                        # discarding triangles out of surface
#                        if self.surface is not None:        
#                            try:
#                                triang.set_mask(~cropSurface(triang, self.surface[:,0], self.surface[:,1]))
#                            except Exception as e:
#                                print('Error in Mesh.show for contouring: ', e)
#                        
#                        self.cax = ax.tricontourf(triang, z, cmap=alpha_map)
#                                                    
                except Exception as e:
                    print('Error in the sensitivity overlay:', e)
        
        if electrodes: #try add electrodes to figure if we have them 
            try: 
#                ax.plot(elec_x, self.elec_z[~iremote],'ko')
                x = np.c_[self.elec_x, self.elec_y, self.elec_z]
                if self.iremote is not None: # it's None for quad mesh
                    x = x[~self.iremote, :]
#                x1 = np.repeat(x, len(x), axis=0)
#                x2 = np.tile(x.T, len(x)).T
#                dist = np.sqrt(np.sum((x1-x2)**2, axis=1))
#                radius = np.nanmin(dist[dist != 0])/3
#                circles = [plt.Circle((xi,yi), radius=radius) for xi,yi in zip(elec_x, self.elec_z[~iremote])]
#                ax.add_collection(PatchCollection(circles, color='k'))
                ax.plot(x[:,0], x[:,2], 'ko', markersize=4)
            except AttributeError:
                print("no electrodes in mesh object to plot")

        # adding interactive display when mouse-over
        centroids = np.array([self.elm_centre[0], self.elm_centre[2]]).T
        def format_coord(x, y):
            dist = np.sqrt(np.sum((centroids - np.array([x, y]))**2, axis=1))
            imin = np.argmin(dist)
            # TODO if imin is out of surface, then don't show value
            return ('x={:.2f} m, elevation={:.2f} m, value={:.3f}'.format(x,y,X[imin]))
        ax.format_coord = format_coord

        # print('Mesh plotted in %6.2f seconds'%(time.time()-t0))
        
        if iplot == True:
            return fig
    
    def draw(self, 
             attr=None,
             edge_color = 'k',
             color_map = 'Spectral',
             color_bar = False,
             vmin= None, vmax = None):
        """
        Redraws a mesh over a prexisting figure canvas, this is intended for saving
        time when plotting each iteration of the resistivity inversion.
        
        Parameters
        ----------
        color_map : string, optional
            color map reference 
        color_bar : boolean, optional 
            `True` to plot colorbar 
        ax : matplotlib axis handle, optional
            Axis handle if preexisting (error will thrown up if not) figure is to be cast to.
        edge_color : string, optional
            Color of the cell edges, set to `None` if you dont want an edge.
        vmin : float, optional
            Minimum limit for the color bar scale.
        vmax : float, optional
            Maximum limit for the color bar scale.
        attr : string, optional
            Which attribute in the mesh to plot, references a dictionary of attributes. `attr` is passed 
            as the key for this dictionary.

        Returns
        ----------
        figure : matplotlib figure 
            Figure handle for the plotted mesh object.
        
        """

        if self.cax == None :
            raise Exception ("mesh canvas variable has not been assigned! Use mesh.show() first")
            
        if attr not in self.attr_cache.keys():
            attr == None
        if attr == None:
            attr = list(self.attr_cache.keys())[3]
            X = np.array(self.attr_cache[attr])
            color_bar_title = attr
        else:
            X = np.array(self.attr_cache[attr])
            color_bar_title = attr
                
        a = time.time() #start timer on how long it takes to plot the mesh
        
        if edge_color == None or edge_color=='none' or edge_color=='None':
            edge_color='face'#set the edge colours to the colours of the polygon patches
            
        if vmin is None:
            vmin = np.min(X)
        if vmax is None:
            vmax = np.max(X)
        
        if attr == 'region':
            cm = plt.get_cmap(color_map, len(np.unique(X)))
            val = np.sort(np.unique(X)).astype(int)
            if len(val) > 1:
                interval = (val[-1]-val[0])/len(val)
                self.cbar.set_ticks(np.arange(val[0]+interval/2, val[-1], interval))
            else:
                self.cbar.set_ticks([1])
            self.cbar.set_ticklabels(val) 
        else:
            cm = color_map
        self.cax.set_cmap(cm) # change the color map if the user wants to 
        
        #following block of code redraws figure 
        self.cax.set_array(X) # set the array of the polygon collection to the new attribute 
        self.cax.set_clim(vmin=vmin, vmax=vmax) # reset the maximum limits of the color map 
        self.ax.add_collection(self.cax)#blit polygons to axis
        self.cbar.set_label(color_bar_title) # change the color bar title 
        self.fig.canvas.draw() # redraw figure canvas (does not make a new figure it is faster fig.show())

        if color_bar:#add the color bar 
           print("you should have decided you wanted a color bar when using the mesh.show function")
            
        print('Mesh plotted in %6.5f seconds'%(time.time()-a))    
    
    
    def show3D(self,color_map = 'Spectral',#displays the mesh using matplotlib
                color_bar = True,
                xlim = None,
                ylim = None,
                zlim = None, 
                ax = None,
                electrodes = True,
                sens = False,
                edge_color = 'k',
                alpha = 1,
                vmax=None,
                vmin=None,
                attr=None,
                elec_color='k',
                use_pyvista=True,
                background_color=(0.8,0.8,0.8),
                pvslices=([],[],[]),
                pvthreshold=None,
                pvgrid=True,
                pvcontour=[]):
        """
        Shows a 3D tetrahedral mesh. 
        
        Parameters
        ----------
        color_map : string, optional
            Matplotlib color map reference 
        color_bar : Boolean, optional 
            `True` to plot colorbar 
        xlim : tuple, optional
            Axis x limits as `(xmin, xmax)`.
        ylim : tuple, optional
            Axis y limits as `(ymin, ymax)`. 
        zlim : tuple, optional
            Axis z limits as `(ymin, ymax)`. 
        ax : matplotlib axis handle, pvista plotter handle, optional
            Axis handle if preexisting (error will thrown up if not) figure is to be cast to.
            If using pyvista then then ax is the plotter the object. 
        electrodes : boolean, optional
            Enter true to add electrodes to plot (if available in mesh class)
        sens : boolean, optional
            Enter true to plot sensitivities (if available in mesh class). Note that for 3D this doesnt work so well. 
        edge_color : string, optional
            Color of the cell edges, set to `None` if you dont want an edge.
        alpha: float, optional
            Should be set between 0 and 1. Sets a transparancy value for the element faces. 
        vmin : float, optional
            Minimum limit for the color bar scale.
        vmax : float, optional
            Maximum limit for the color bar scale.
        attr : string, optional
            Which attribute in the mesh to plot, references a dictionary of attributes. attr is passed 
            as the key for this dictionary.
        elec_color : string, optional
            Colour of the electrodes on the plot if electrodes = True. Default
            is 'k' for black. Can be a 3 by 1 tuple or string identifier. 
        use_pyvista : bool, optional
            Use visual toolkit backend for displaying 3D mesh, note that pyvista
            must be installed for this to work. 
        background_color : tuple, optional 
            Background color assigned to pyvista plotter object when created. Not yet
            supported for matplotlib axis handles. 
        pvslices : tuple of list of float, optional
            Determine the X, Y, Z slices. e.g.: ([3], [], [-3, -4]) will add
            a slice normal to X in 3 and two slices normal to Z in -3 and -4.
        pvthreshold : list of two floats, optional
            Keep values between pvthreshold[0] and pvthreshold[1].
        pvgrid : bool, optional
            Show grid or not.
        pvcontour: list of float, optional
            Values of the isosurface to be plotted.

        Returns
        -------
        figure : matplotlib figure
            Figure handle for the plotted mesh object.
        
        Notes
        -----
        Show a mesh object using matplotlib. The color map variable should be 
        a string refering to the color map you want (default is "Spectral").
        As we're using the matplotlib package here any color map avialable within 
        matplotlib package can be used to display the mesh here also. See: 
        https://matplotlib.org/2.0.2/examples/color/colormaps_reference.html
        
        Plotting sensitivies using sens=True is not reccomended. The matplotlib renderer has trouble with it. 
        """
        if not isinstance(color_map,str):#check the color map variable is a string
            raise NameError('color_map variable is not a string')
        
        #decide which attribute to plot, we may decide to have other attritbutes! 
        if attr not in self.attr_cache.keys():
            attr == None
        if attr == None:
            attr = list(self.attr_cache.keys())[3]
            X = np.array(self.attr_cache[attr])
            color_bar_title = attr
        else:
            X = np.array(self.attr_cache[attr])
            color_bar_title = attr
                
        #determine how to crop the mesh  
        if zlim == None:
            zlim = [min(self.node_z), max(self.node_z)]
        if self.iremote is not None: 
            if xlim == None:
                xlim = [min(self.elec_x[~self.iremote]), max(self.elec_x[~self.iremote])]
            if ylim == None:
                ylim = [min(self.elec_y[~self.iremote]), max(self.elec_y[~self.iremote])]
        else:
            try: 
                if xlim == None:
                    xlim = [min(self.elec_x), max(self.elec_x)]
                if ylim == None:
                    ylim = [min(self.elec_y), max(self.elec_y)]
            except AttributeError: #if no electrodes present use the node limits 
                if xlim == None:
                    xlim = [min(self.node_x), max(self.node_x)]
                if ylim == None:
                    ylim = [min(self.node_y), max(self.node_y)]
           
        if abs(xlim[0] - xlim[1]) < 0.001:# protection against thin axis margins 
            xlim=[xlim[0]-2,xlim[1]+2]
        if abs(zlim[0] - zlim[1]) < 0.001:
            zlim=[zlim[0]-2,zlim[1]+2]
        if abs(ylim[0] - ylim[1]) < 0.001:
            ylim=[ylim[0]-2,ylim[1]+2]
            
        if vmin is None: # determine min and max of X array 
            vmin = np.min(X)
        if vmax is None:
            vmax = np.max(X)
        
        #### use pyvista for 3D viewing ####         
        if use_pyvista and pyvista_installed:
            nmesh = self.copy() # new mesh which is cropped down 
            nmesh.attr_cache = {color_bar_title:X} # make the attr of interest the only attribute
            folder = tempfile.TemporaryDirectory()
            fname = os.path.join(folder.name, '__to_pv_mesh.vtk')
            nmesh.write_vtk(fname)
            self.pvmesh = pv.read(fname)
            folder.cleanup()
            #### possible to iniate pyvista mesh directly with vertices and connection matrix #### 
            #### however the below code is apparently unstable on pyvista 0.24.0 #### 
            # self.pvmesh = pv.PolyData(np.array([self.node_x,self.node_y,self.node_z]).T, 
            #                           np.array(self.con_matrix).T)
            # self.pvmesh[color_bar_title] = X
            
            #crop down to the bounding box 
            self.pvmesh = self.pvmesh.clip_box((xlim[0],xlim[1],ylim[0],ylim[1],zlim[0],zlim[1]),invert=False)
                        
            if edge_color is None or edge_color=='none' or edge_color=='None':
                edges = False # then dont show element edges 
            else:
                edges = True
            if ax is None: # make a plotter object if not already given 
                ax = pv.BackgroundPlotter()
                ax.background_color = background_color
            else: # check the ax argument is for pyvista not matplotlib 
                typ_str = str(type(ax))
                if typ_str.find('pyvista') == -1:
                    raise Exception('Error plotting with pyvista, show3D (meshTools.py) expected a pyvista plotter object but got %s instead'%typ_str)
                ax.set_background(background_color)
                
            if pvthreshold is not None:
                if isinstance(pvthreshold, list):
                    if pvthreshold[0] is None:
                        pvthreshold[0] = np.nanmin(X)
                    if pvthreshold[1] is None:
                        pvthreshold[1] = np.nanmax(X)
                self.pvmesh = self.pvmesh.threshold(value=pvthreshold)
            if len(pvcontour) > 0:
                self.pvmesh = self.pvmesh.cell_data_to_point_data()
                self.pvmesh = self.pvmesh.contour(isosurfaces=pvcontour)
            if pvgrid:
                ax.show_grid(color='k')
                
            if np.sum([len(a) for a in pvslices]) > 0: # we have slices
                ax.add_mesh(self.pvmesh.outline(), color='k')
                for i, ss in enumerate(pvslices): # X, Y then Z slice
                    normal = np.zeros(3)
                    normal[i] = 1
                    for s in ss:
                        if ((s > np.nanmin(self.elm_centre[i])) &
                            (s < np.nanmax(self.elm_centre[i]))):
                            origin = np.zeros(3)
                            origin[i] = s
                            mesh_slice = self.pvmesh.slice(normal=normal, origin=origin)
                            if mesh_slice.number_of_points > 0 or mesh_slice.number_of_cells > 0:
                                ax.add_mesh(mesh_slice,
                                            cmap=color_map,
                                            clim=[vmin, vmax],
                                            show_scalar_bar=color_bar,
                                            show_edges=edges,
                                            opacity=alpha,
                                            scalar_bar_args={'color':'k'})
                            else:
                                print('empty mesh')
            else:        
                if self.pvmesh.number_of_points > 0 or self.pvmesh.number_of_cells > 0:
                    ax.add_mesh(self.pvmesh,
                                cmap=color_map, #matplotlib colormap 
                                clim=[vmin,vmax], #color bar limits 
                                show_scalar_bar=color_bar,#plot the color bar? 
                                show_edges=edges, #show edges
                                opacity=alpha,
                                scalar_bar_args={'color':'k',# 'interactive':True,
                                                 'vertical':False,
                                                 'title_font_size':16,
                                                 'label_font_size':14})
                else:
                    print('empty mesh')
                    
            if electrodes: # then add the electrodes to the plot 
                try:
                    points = np.array([self.elec_x,self.elec_y,self.elec_z]).T
                    pvelec = pv.PolyData(points)
                    ax.add_mesh(pvelec, color=elec_color, point_size=10.,
                                render_points_as_spheres=True)
                except AttributeError as e:
                    print("Could not plot 3d electrodes, error = "+str(e))
            
            ax.show()
            return
        #### else fall back on to matplotlib scheme ####    
        
        #if prism mesh use that function instead 
        if self.ndims==2:
            warnings.warn("Its reccomended to use mesh.show() for 2D meshes, results of 3D show will be unstable")
        elif self.type2VertsNo() == 6: # use column mesh show instead 
            self.show_prism_mesh()
            return 
        
        t0 = time.time() # benchmark function
        
        #make 3D figure 
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        
        #set axis limits     
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_zlim(zlim) # doesn't seem to work in the UI
        #set color bar limits
            
        if edge_color == None or edge_color=='none' or edge_color=='None':
            edge_color='face'#set the edge colours to the colours of the polygon patches
        
        # search through each element to see if it is on the edge of the mesh, 
        # this step is important as it is very expensive to plot anything in 3D using matplotlib 
        # triangles on the edge of the mesh will be used only once
        
        elm_x = self.elm_centre[0]
        elm_y = self.elm_centre[1]
        elm_z = self.elm_centre[2]
        in_elem = in_box(elm_x,elm_y,elm_z,xlim[1],xlim[0],ylim[1],ylim[0],zlim[1],zlim[0])#find elements veiwable in axis
        X = X[in_elem]#reassign X to elements inside the box limits 
        temp_con_mat = np.array(self.con_matrix,dtype='int64')#temporary connection matrix which is just the elements inside the box
        con_mat=temp_con_mat[:,in_elem] # truncate elements
        inside_numel = len(con_mat[0])#number of inside elements 
        tri_combo = np.zeros((inside_numel,4),dtype='float64')
        
        S = [0]*self.num_elms 
        if sens:
            try:
                S = self.sensitivities
            except AttributeError:
                sens = False
                print('no sensitivities to plot')
        
        for i in range(inside_numel):
            idx1 = con_mat[0][i]#extract indexes 
            idx2 = con_mat[1][i]
            idx3 = con_mat[2][i]
            idx4 = con_mat[3][i]
            
            face1 = idx1*idx2*idx3 # assign each face a code
            face2 = idx1*idx2*idx4
            face3 = idx2*idx3*idx4
            face4 = idx1*idx4*idx3

            tri_combo[i,0] = face1#face 1 
            tri_combo[i,1] = face2#face 2 
            tri_combo[i,2] = face3#face 3 
            tri_combo[i,3] = face4#face 4 
            
        temp,index,counts = np.unique(tri_combo.flatten(),return_index=True,return_counts=True) # find the unique values 
        single_vals_idx = counts==1 # faces on edge of volume only appear once
        edge_element_idx = index[single_vals_idx]/4
        face_element_idx = np.floor(edge_element_idx)
        face_probe = edge_element_idx - np.floor(edge_element_idx)
        
        truncated_numel = len(face_element_idx)
        face_list = [0] * truncated_numel
        assign = [0] * truncated_numel # the number assigned to each face
        sensi = [0] * truncated_numel # the sensitivity assigned to each face
        
        for i in range(truncated_numel):
            ref = int(face_element_idx[i])
            idx1 = con_mat[0][ref]
            idx2 = con_mat[1][ref]
            idx3 = con_mat[2][ref]
            idx4 = con_mat[3][ref]
            
            vert1 = (self.node_x[idx1],self.node_y[idx1],self.node_z[idx1])
            vert2 = (self.node_x[idx2],self.node_y[idx2],self.node_z[idx2])
            vert3 = (self.node_x[idx3],self.node_y[idx3],self.node_z[idx3])
            vert4 = (self.node_x[idx4],self.node_y[idx4],self.node_z[idx4])
            
            if face_probe[i] == 0: #if single_val_idx. == 0 > face1
                face_list[i] = (vert1,vert2,vert3)#face 1 
            elif face_probe[i] == 0.25:#if single_val_idx. == 0.25 > face2
                face_list[i] = (vert1,vert2,vert4)#face 2
            elif face_probe[i] == 0.5:#if single_val_idx. == 0.5 > face3
                face_list[i] = (vert2,vert3,vert4)#face 3 
            elif face_probe[i] == 0.75:#if single_val_idx. == 0.75 > face4
                face_list[i] = (vert1,vert4,vert3)#face 4
            
            assign[i] = X[ref]#get attribute value into assigned array
            sensi[i] = S[ref]
          
        polly = Poly3DCollection(face_list,linewidth=0.5) # make 3D polygon collection
        polly.set_alpha(alpha)#add some transparancy to the elements
        try:
            polly.set_array(np.array(assign))
        except MemoryError:#catch this error and print something more helpful than matplotlibs output
            raise MemoryError("Memory access voilation encountered when trying to plot mesh, \n please consider truncating the mesh or display the mesh using paraview.")
        try:
            polly.set_edgecolor(edge_color) #seems to be a bug with matplotlib where this function doesnt work in 3.1
        except AttributeError:
            pass 
        
        polly.set_cmap(color_map) # set color map 
        polly.set_clim(vmin=vmin, vmax=vmax) # reset the maximum limits of the color map 
        ax.add_collection3d(polly, zs=0, zdir='z') # for matplotlib > 3.0.2
        self.cax = polly
        
        if color_bar:#add the color bar 
            self.cbar = plt.colorbar(self.cax, ax=ax, format='%.1f')
            self.cbar.set_label(color_bar_title) #set colorbar title
            
#        ax.set_aspect('equal')#set aspect ratio equal (stops a funny looking mesh)
        
        if sens: #add sensitivity to plot if available
            weights = np.array(sensi) #values assigned to alpha channels 
            alphas = np.linspace(1, 0, len(sensi))#array of alpha values 
            raw_alpha = np.ones((len(sensi),4),dtype=float) #raw alpha values 
            raw_alpha[..., -1] = alphas
            alpha_map = ListedColormap(raw_alpha) # make a alpha color map which can be called by matplotlib
            #make alpha collection
            alpha_coll = Poly3DCollection(face_list, array=weights, cmap=alpha_map, edgecolors='none', linewidths=0)#'face')
            ax.add_collection(alpha_coll)
            
        if electrodes: #try add electrodes to figure if we have them 
            try: 
                ax.scatter(self.elec_x,self.elec_y,zs=np.array(self.elec_z),
                           s=20, c=elec_color, marker='o')
                #note 3D plotting doesnt work well in matplotlib, use of 
                #pyvista is recomended here. 
            except AttributeError as e: # no electrodes in mesh object. 
                print("Could not plot 3d electrodes, error = "+str(e))
            
        print('Mesh plotted in %6.5f seconds'%(time.time()-t0))



    def showSlice(self, attr='Resistivity(log10)', axis='z', vmin=None, vmax=None, ax=None):
        """ Show 3D mesh slice.
        
        Parameters
        ----------
        attr : str, optional
            String the of the variable to plot. Default is `Resistivity(log10)`.
        axis : str, optional
            Axis to be sliced either, x, y or z.
        vmin : float, optional
            Minimum value for colorbar.
        vmax : float, optional
            Maximum value for colorbar.
        ax : matplotlib.Axis
            If provided, plot will be drawn on this axis.
        """
        values = np.array(self.attr_cache[attr])        
        dimDico = {'x':0,'y':1,'z':2}
        dim = dimDico[axis]
        elms = np.array(self.con_matrix).T
        nodes = np.array([self.node_x, self.node_y, self.node_z]).T
        sliceMesh(nodes, elms, values, label=attr, dim=dim, vmin=vmin, vmax=vmax, ax=ax)
        
        
    def show_prism_mesh(self,color_map = 'Spectral',#displays the mesh using matplotlib
             color_bar = True,
             ax = None,
             electrodes = True,
             sens = False,
             edge_color = 'k',
             alpha = 1,
             aspect = 'equal',
             vmax=None,
             vmin=None,
             attr=None,
             xlim = None,
             ylim = None,
             zlim = None):
        """
        Shows a 3D prism mesh. 
        
        Parameters
        ----------
        color_map : string, optional
            color map reference 
        color_bar : Boolean, optional 
            `True` to plot colorbar 
        xlim : tuple, optional
            Axis x limits as `(xmin, xmax)`.
        ylim : tuple, optional
            Axis y limits as `(ymin, ymax)`. 
        zlim : tuple, optional
            Axis z limits as `(ymin, ymax)`. 
        ax : matplotlib axis handle, optional
            Axis handle if preexisting (error will thrown up if not) figure is to be cast to.
        electrodes : boolean, optional
            Enter true to add electrodes to plot (if available in mesh class)
        sens : boolean, optional
            Enter true to plot sensitivities (if available in mesh class). Note that for 3D this doesnt work so well. 
        edge_color : string, optional
            Color of the cell edges, set to `None` if you dont want an edge.
        alpha: float, optional
            Should be set between 0 and 1. Sets a transparancy value for the element faces. 
        vmin : float, optional
            Minimum limit for the color bar scale.
        vmax : float, optional
            Maximum limit for the color bar scale.
        attr : string, optional
            Which attribute in the mesh to plot, references a dictionary of attributes. attr is passed 
            as the key for this dictionary.
        aspect = string, optional
            defines the aspect ratio of the plot.
            'equal' locks the aspect ratio.
            'auto', aspect ratio is define by plotting area.

        Returns
        ----------
        figure : matplotlib figure 
            Figure handle for the plotted mesh object.
        """
        
        if not isinstance(color_map,str):#check the color map variable is a string
            raise NameError('color_map variable is not a string')
        
        if self.ndims==2:
            warnings.warn("Use Mesh.show() for 2D meshes")
        elif self.type2VertsNo() != 6:
            warnings.warn("Function can only be used to display prism meshes")
            return 
            
        #decide which attribute to plot, we may decide to have other attritbutes! 
        if attr not in self.attr_cache.keys():
            attr == None
        if attr == None:
            attr = list(self.attr_cache.keys())[3]
            X = np.array(self.attr_cache[attr])
            color_bar_title = attr
        else:
            X = np.array(self.attr_cache[attr])
            color_bar_title = attr

        t0 = time.time() # benchmark function
        
        #make 3D figure 
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        
        if zlim==None:
            zlim=[min(self.node_z),max(self.node_z)]
        if xlim==None:
            xlim=[min(self.node_x), max(self.node_x)]
        if ylim==None:
            ylim=[min(self.node_y), max(self.node_y)]
        #set axis limits     
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_zlim(zlim) # doesn't seem to work in the UI
        #set color bar limits

        if vmin is None:
            vmin = np.min(X)
        if vmax is None:
            vmax = np.max(X)
            
        if edge_color == None or edge_color=='none' or edge_color=='None':
            edge_color='face'#set the edge colours to the colours of the polygon patches
            
        #construct patches 
        #print('constructing patches')
        S = [0]*self.num_elms 
        if sens:
            try:
                S = self.sensitivities
            except AttributeError:
                sens = False
                print('no sensitivities to plot')
                
        con_mat = self.con_matrix # just plotting top and bottom parts of elements 
        combo = np.zeros((self.num_elms,5),dtype='float64')
        
        for i in range(self.num_elms):
            idx1 = con_mat[0][i]
            idx2 = con_mat[1][i]
            idx3 = con_mat[2][i]
            idx4 = con_mat[3][i]
            idx5 = con_mat[4][i]
            idx6 = con_mat[5][i]
            
            #assign each face a code 
            fc1 = int(str(idx1)+str(idx2)+str(idx3))#face code 1 
            fc2 = int(str(idx4)+str(idx5)+str(idx6))#face code 2 
            fc3 = int(str(idx1)+str(idx2)+str(idx5)+str(idx4))#face code 3
            fc4 = int(str(idx2)+str(idx3)+str(idx6)+str(idx5))#face code 4 
            fc5 = int(str(idx1)+str(idx3)+str(idx6)+str(idx4))#face code 5
            
            combo[i,0] = fc1
            combo[i,1] = fc2
            combo[i,2] = fc3
            combo[i,3] = fc4
            combo[i,4] = fc5
            
        #return combo
        temp,index,counts = np.unique(combo.flatten(),return_index=True,return_counts=True) 
        single_vals_idx = counts==1 # faces on edge of volume only appear once
        edge_element_idx = index[single_vals_idx]/5
        face_element_idx = np.floor(edge_element_idx)
        face_probe = edge_element_idx - np.floor(edge_element_idx)
        face_probe = np.round(face_probe,1) # round to nearest 1 decimal 
        
        
        truncated_numel = len(face_element_idx)
        face_list = [()] * truncated_numel
        assign = [0] * truncated_numel # the number assigned to each face
        sensi = [0] * truncated_numel # the sensitivity assigned to each face
            
        for i in range(truncated_numel):
            ref = int(face_element_idx[i])
            idx1 = con_mat[0][ref]
            idx2 = con_mat[1][ref]
            idx3 = con_mat[2][ref]
            idx4 = con_mat[3][ref]
            idx5 = con_mat[4][ref]
            idx6 = con_mat[5][ref]
            vert1 = (self.node_x[idx1],self.node_y[idx1],self.node_z[idx1])
            vert2 = (self.node_x[idx2],self.node_y[idx2],self.node_z[idx2])
            vert3 = (self.node_x[idx3],self.node_y[idx3],self.node_z[idx3])
            vert4 = (self.node_x[idx4],self.node_y[idx4],self.node_z[idx4])
            vert5 = (self.node_x[idx5],self.node_y[idx5],self.node_z[idx5])
            vert6 = (self.node_x[idx6],self.node_y[idx6],self.node_z[idx6])
            
            if face_probe[i] == 0: #if single_val_idx. == 0 > face1
                face_list[i] = (vert1,vert2,vert3)#face 1 
            elif face_probe[i] == 0.2:#if single_val_idx. == 0.2 > face2
                face_list[i] = (vert4,vert5,vert6)#face 2
            elif face_probe[i] == 0.4:#if single_val_idx. == 0.4 > face3
                face_list[i] = (vert1,vert2,vert5,vert4)#face 3 
            elif face_probe[i] == 0.6:#if single_val_idx. == 0.6 > face4
                face_list[i] = (vert2,vert3,vert6,vert5)#face 4
            elif face_probe[i] == 0.8:#if single_val_idx. == 0.8 > face5
                face_list[i] = (vert1,vert3,vert6,vert4)
                
            sensi[i] = S[ref]
            assign[i] = X[ref]
            
        #add patches to 3D figure 
        polly = Poly3DCollection(face_list,linewidth=0.5) # make 3D polygon collection
        polly.set_alpha(alpha)#add some transparancy to the elements
        try:
            polly.set_array(np.array(assign))
        except MemoryError:#catch this error and print something more helpful than matplotlibs output
            raise MemoryError("Memory access voilation encountered when trying to plot mesh, \n please consider truncating the mesh or display the mesh using paraview.")
        polly.set_edgecolor(edge_color)
        polly.set_cmap(color_map) # set color map 
        polly.set_clim(vmin=vmin, vmax=vmax) # reset the maximum limits of the color map 
#        ax.add_collection3d(polly, zs='z')#blit polygons to axis 
        ax.add_collection3d(polly, zs=0, zdir='z') # for matplotlib > 3.0.2
        self.cax = polly
        
        if color_bar:#add the color bar 
            self.cbar = plt.colorbar(self.cax, ax=ax, format='%.1f')
            self.cbar.set_label(color_bar_title) #set colorbar title
            
        # ax.set_aspect(aspect)#set aspect ratio equal (stops a funny looking mesh)
        
        if sens: #add sensitivity to plot if available
            weights = np.array(sensi) #values assigned to alpha channels 
            alphas = np.linspace(1, 0, len(sensi))#array of alpha values 
            raw_alpha = np.ones((len(sensi),4),dtype=float) #raw alpha values 
            raw_alpha[..., -1] = alphas
            alpha_map = ListedColormap(raw_alpha) # make a alpha color map which can be called by matplotlib
            #make alpha collection
            alpha_coll = Poly3DCollection(face_list, array=weights, cmap=alpha_map, edgecolors='none', linewidths=0)#'face')
            ax.add_collection(alpha_coll)
            
        if electrodes: #try add electrodes to figure if we have them 
            try: 
                ax.scatter(self.elec_x,self.elec_y,zs=np.array(self.elec_z),
                           s=20, c='k', marker='o')
            except AttributeError as e:
                print("could not plot 3d electrodes, error = "+str(e))
                
        print('Mesh plotted in %6.5f seconds'%(time.time()-t0))
        
    #%% mesh zoning     
    def assign_zone(self,poly_data):
        """ Assign material/region assocations with certain elements in the mesh 
        say if you have an area you'd like to forward model. 
        ***2D ONLY***
            
        Parameters
        ----------
        poly_data : list 
            list of 2 by N (x,y) numpy arrays describing each polygon.
            
        Returns
        ---------
        material_no : numpy.array
            Element associations starting at 1. So 1 for the first region 
            defined in the region_data variable, 2 for the second region 
            defined and so on. If the element can't be assigned to a region
            then it'll be left at 0. 
        """   
        if self.ndims==3:
            raise ValueError('Assign zone available with 2D meshes only')
            
        no_elms=self.num_elms#number of elements 
        elm_xz=np.array([self.elm_centre[0],self.elm_centre[2]]).T#centriods of 2d mesh elements 
        material_no=np.zeros(no_elms,dtype=int)#attribute number
        
        if not isinstance(poly_data,dict):
            raise Exception("poly_data input is not a dictionary")
        
        #now on to extracting the data of interest
        print('Assigning element attribute IDs...')
        for i in range(len(poly_data)):
            polyline = poly_data[i]
            path = mpath.Path(polyline)
            inside = path.contains_points(elm_xz) 
            material_no[inside]=i+1
                            
        self.zone = material_no
        self.add_attribute(material_no,'zone')
        return material_no     
        
    def assign_zone_attribute(self,material_no,attr_list,new_key):
        """ Asssigns values to the mesh which depend on region / material only. E.G 
        a single resistivity value.
            
        Parameters
        ----------
        material_no : array or list
            Integers starting at 0 or 1, and ascend in intervals of 1, which 
            correspond to a material in the mesh returned from assign_attr_ID.
            Should have the same length as the number of elements in the mesh.
        attr_list : list
            Values corresponding to a material number in the mesh. eg. if you had 3 regions in the mesh then you give
            `[resistivity1,resistivity2,resistivity3]`.
        new_key : string
            Key identifier assigned to the attribute in the attr_cache. 
        
        Notes  
        -----
        Mesh object will now have the new attribute added once the function is run.
        Use the `mesh.show()` (or `.draw()`) function to see the result. 
        """ 
        if len(material_no) != self.num_elms:
            raise ValueError("Mismatch between the number of elements and material propeties")
        
        new_para=np.array([0]*self.num_elms)
        
        if min(material_no)==1:#cor_fac allows for compatability with an index system starting at 1 or 0 
            cor_fac=1
        else:
            cor_fac=0
            
        material_no = np.array(material_no) - cor_fac
        
        for i in range(len(attr_list)):
            idx = material_no == i
            new_para[idx] = attr_list[i]
        
        self.attr_cache[new_key] = new_para
        self.no_attributes += 1            

    def apply_func(self,mesh_paras,material_no,new_key,function,*args):
        """ Applies a function to a mesh by zone number and mesh parameter.
        
        Parameters
        ----------
        mesh_paras : array like
            Mesh parameters from which new parameters are calculated.
        material_no : array like of ints
            Material type assigned to each element, should be numbered consectively from 1 to n. in the form 1 : 1 : 2 : n.
            ...ie if you have 2 materials in the mesh then pass an array of ones and twos. zeros will be ignored. 
        new_key : string
            Key assigned to the parameter in the attr_cache. DOES NOT default.
        function : function
            Function to be applied to mesh attributes, first argument must be the mesh parameter.
        args : [see function info]
            All arguments to be passed through function after to modify the mesh parameter,
            ... argument must be in the form of [(argA1,argB1),(argA2,argB2)], 
            ... where letters are the material, numbers refer to the argument number
        
        Notes  
        ---------
        Mesh object will now have the new attribute added once the function is run.
        Use the `mesh.show()` (or `.draw()`) function to see the result. 
        """
    
        if len(material_no)!=len(mesh_paras):
            raise ValueError('Mismatch between the number of material propeties (for the mesh) and parameters to be converted')
        new_para=[0]*self.num_elms
        #iterate through each set of argument variables
        for iteration in range(len(args[0])):
            parameters=[items[iteration] for items in args]#return parameters 
            parameters.insert(0,0)#this adds an element to the front of the parameters which can be swapped out to mesh_paras
            for i in range(self.num_elms):
                if material_no[i]==iteration+1:#does the material match the iteration? 
                    parameters[0]=mesh_paras[i]#change parameter value at start of variables list
                    new_para[i]=function(*parameters)#compute new parameter   
        self.attr_cache[new_key] = new_para
        self.no_attributes += 1
        #return new_para
        
    def computeElmDepth(self):
        """ Compute the depth of elements relative to the surface of the mesh.
        """
        if self.ndims == 2: # use 1D interpolation
            xz = self.extractSurface()
            datum_x = xz[0]
            datum_z = xz[1]
            
            elm_x = np.array(self.elm_centre[0])
            elm_z = np.array(self.elm_centre[2])
            min_idx = np.argmin(datum_x)
            max_idx = np.argmax(datum_x)
            Z = np.interp(elm_x,datum_x,datum_z,left=datum_z[min_idx],
                          right=datum_z[max_idx])
            depth = Z - elm_z
            self.attr_cache['depths'] = depth
            self.no_attributes += 1
            return depth
        if self.ndims == 3: # use 2D interpolation
            elm_x = np.array(self.elm_centre[0])
            elm_y = np.array(self.elm_centre[1])
            elm_z = np.array(self.elm_centre[2])  
            #use interpolation to work out depth to datum 
            surmesh = self.extractSurface()
            Z = interp.triangulate(elm_x, elm_y, surmesh.node_x, 
                                   surmesh.node_y, surmesh.node_z)
            depth = Z - elm_z
            self.attr_cache['depths'] = depth # add cell depths to attribute cache
            self.no_attributes += 1
            return depth
            
    def moveElecNodes(self, new_x, new_y, new_z, debug=True):
        """ Move the electrodes to different nodes which are close to the given coordinates. 
        This is useful for timelapse surveys where the electrodes move through time, 
        ideally this is implimented on a mesh which is refined near the surface. If 
        no nodes are assigned to mesh, a mesh.e_nodes variable is created.  
        
        Parameters
        ------------
        new_x: array l ike
            new electrode x coordinates 
        new_y: array like
            new electrode y coordinates, if 2D mesh let new_y=None and the array will automatically be
            assigned an array of zeros. 
        new_z: array-like
            new electrode z coordinates 
        debug: bool, optional
            Controls if any changes to electrode nodes will be output to console. 
        Returns
        ------------
        node_in_mesh: np array
            Array of mesh node numbers corresponding to the electrode postions/
            coordinates.
            
        Notes
        ------------
        Mesh.e_nodes parameter is updated (or added) after running function. 
        """
        #formalities 
        if len(new_x) != len(new_y) and len(new_x) != len(new_z):#set up protection
            raise ValueError("mismatch in the length of the new_x, new_y and new_z arrays")
        try:
            if len(new_x) != len(self.elec_x) or len(new_z) != len(self.e_nodes):
                raise ValueError("mismatch between the length of new electrode position array and old array")
            has_nodes = True
        except AttributeError:
            has_nodes = False
            
        if new_y is None:
            new_y = np.zeros_like(new_x)
            
        self.node_x = np.array(self.node_x)
        self.node_y = np.array(self.node_y)
        self.node_z = np.array(self.node_z)
            
        node_in_mesh = [0]*len(new_x)    
        for i in range(len(new_x)):
            sq_dist = (self.node_x - new_x[i])**2 + (self.node_y - new_y[i])**2 + (self.node_z - new_z[i])**2 # find the minimum square distance
            node_in_mesh[i] = np.argmin(sq_dist) # min distance should be zero, ie. the node index.
            if has_nodes and debug:
                if node_in_mesh[i] != self.e_nodes[i]:
                    print("Electrode %i moved from node %i to node %i"%(i,node_in_mesh[i],self.e_nodes[i]))#print to show something happening
        
        self.add_e_nodes(node_in_mesh) # update e_node parameter
        if len(np.unique(node_in_mesh)) != len(new_x):
            warnings.warn("The number of new electrode nodes does not match the number of electrodes, which means a duplicated node is present! Please make mesh finer.")   
     
    
        return np.array(node_in_mesh, dtype=int) # note this is the node position with indexing starting at 0. 


    def write_dat(self, file_path='mesh.dat'):
        """Write a mesh.dat kind of file for mesh input for R2. R2 takes a mesh
        input file for triangular and quadrilateral meshes.
        
        Parameters
        ----------
        file_path : str, optional
            Path to the file. By default 'mesh.dat' is saved in the working directory.
        """
        if not isinstance(file_path,str):
            raise TypeError("expected string argument for file_path")
        ### write data to mesh.dat kind of file ###
        #open mesh.dat for input      
        with open(file_path, 'w') as fid:
            
            #write to mesh.dat total num of elements and nodes
            # find furthest node from first electrode to be dirichlet node
            xyz = np.array([self.node_x, self.node_y, self.node_z]).T # N x 3
            if self.e_nodes is not None:
                idirichlet = np.argmax(np.sqrt(np.sum((xyz-xyz[self.e_nodes[0],:])**2, axis=1)))
            else:
                idirichlet = len(self.node_x)
            if self.ndims == 3:
                fid.write('%i %i %i 0 %i\n'%(self.num_elms,self.num_nodes,1,self.type2VertsNo()))
            else:
                fid.write('%i %i %i\n'%(self.num_elms,self.num_nodes,idirichlet))
            param = np.array(self.attr_cache['param'])
            zone = np.array(self.attr_cache['zones'])
    
            #write out elements         
            no_verts = self.type2VertsNo()
            for i in range(self.num_elms):
                elm_no=i+1
                fid.write("%i "%elm_no)
                [fid.write("%i "%(self.con_matrix[k][i]+1)) for k in range(no_verts)]
                fid.write("%i %i\n"%(param[i],zone[i]))
    
            #now add nodes
            if self.ndims == 3:
                for i in range(self.num_nodes):
                    ni_no=i+1
                    fid.write("%i %16.8f %16.8f %16.8f\n"%#node number, x coordinate, y coordinate, z coordinate
                              (ni_no,
                               self.node_x[i],
                               self.node_y[i],
                               self.node_z[i]))
                fid.write('{:d}'.format(idirichlet))
            else:
                for i in range(self.num_nodes):
                    ni_no=i+1
                    fid.write("%i %16.8f %16.8f\n"%#node number, x coordinate, y coordinate
                              (ni_no,
                               self.node_x[i],
                               self.node_z[i]))


    def write_vtk(self, file_path="mesh.vtk", title=None, replace_nan=-9999):
        """Writes a vtk file for the mesh object, everything in the attr_cache
        will be written to file as attributes. We suggest using Paraview 
        to display the mesh outside of ResIPy. It's fast and open source :). 
        
        Parameters
        ----------
        file_path: str, optional
            Maps where python will write the file, if left as `default` then mesh.vtk
            will be written the current working directory. 
        title : str, optional
            Header string written at the top of the vtk file .
        """
        #formalities 
        if title == None:
            try:
                title = self.mesh_title
            except AttributeError:
                title = "output from resipy meshTools module"
        if not file_path.endswith('.vtk'):
            file_path +='.vtk'#append .vtk extension to end of file path if not there
        #open file and write header information  
        fh = open(file_path,'w')
        fh.write("# vtk DataFile Version 3.0\n")
        fh.write(title+"\n")
        fh.write("ASCII\nDATASET UNSTRUCTURED_GRID\n")
        #define node coordinates
        fh.write("POINTS %i double\n"%self.num_nodes)
        for i in range(self.num_nodes):
            fh.write("%8.6f\t%8.6f\t%8.6f\n"%(self.node_x[i],self.node_y[i],self.node_z[i]))
        #define the connection matrix    
        no_verts = self.type2VertsNo()
        no_readable = self.num_elms*(1+no_verts)
        fh.write("CELLS %i %i\n"%(self.num_elms,no_readable))
        for i in range(self.num_elms):
            fh.write("%i\t"%no_verts)
            for k in range(no_verts):
                fh.write("{}    ".format(self.con_matrix[k][i]))
            fh.write("\n")
        #cell types
        fh.write("CELL_TYPES %i\n"%self.num_elms)
        [fh.write("%i "%self.cell_type[0]) for i in range(self.num_elms)];fh.write("\n")
        #write out the data
        fh.write("CELL_DATA %i\n"%self.num_elms)
        attr_cache = self.attr_cache
        for i,key in enumerate(attr_cache.keys()):
            fh.write("SCALARS %s double 1\n"%key.replace(' ','_'))
            fh.write("LOOKUP_TABLE default\n")
            X = np.array(attr_cache[key])
            X[np.isnan(X)]=replace_nan
            [fh.write("%8.6f "%X[j]) for j in range(self.num_elms)]
            fh.write("\n")
        
        #finish writing
        fh.write("POINT_DATA %i"%self.num_nodes)        
        fh.close()
    

    def write_attr(self,attr_key,file_name='_res.dat'):
        """ Writes a attribute to a _res.dat type file. file_name entered
        seperately because it will be needed for the R2 config file.
        The reason for this function is so you can write a forward model 
        parameter input file. 
        
        Parameters
        ----------
        attr_key: string
            Key identifying the attr to be written in the mesh object attr_cache.
        file_name: string, optional
            Name of the _res.dat type file.
        """
        #formality checks
        #double check the fname is not longer than 15 characters as the code doesnt like it apparently 
        if len(ntpath.basename(file_name))>15:
            raise NameError("File name for _res.dat type file cannot be longer than 15 characters")
            
        if isinstance(file_name,str)==False:
            raise NameError("file_name argument must be a string")
        
        #the format of the _res.dat file is such that
        #| x coordinate | y coordinate | value | log(value) | 
        fh = open(file_name,'w')#open file handle 
        x_coords=self.elm_centre[0]#get element coordinates
        y_coords=self.elm_centre[1]
        z_coords=self.elm_centre[2]
        values=self.attr_cache[attr_key]
        log_values=np.log10(np.array(values))
        if self.ndims==3:
            for i in range(self.num_elms):
                fh.write("\t{: 10.5e}\t{: 10.5e}\t{: 10.5e}\t{: 10.5e}\t{: 10.5e}\n".format(x_coords[i],y_coords[i],z_coords[i],values[i],log_values[i]))
        else:
            for i in range(self.num_elms):
                fh.write("\t{: 10.5e}\t{: 10.5e}\t{: 10.5e}\t{: 10.5e}\n".format(x_coords[i],z_coords[i],values[i],log_values[i]))
            
        fh.close()
        
    def write_csv(self,file_name='mesh.csv'):
        """ Write a .csv file of the mesh, the first 3 columns are the element 
        centres at coordinates x,y,z and the rest of the columns are the 
        attributes in the attr_cache

        Parameters
        ----------
        file_name : String, optional
            The default is 'mesh.csv'.

        Returns
        -------
        None.

        """
        if isinstance(file_name,str)==False:
            raise NameError("file_name argument must be a string")
        x_coords=self.elm_centre[0]#get element coordinates
        y_coords=self.elm_centre[1]
        z_coords=self.elm_centre[2]
        keys = self.attr_cache.keys()
        
        #open file 
        fh = open(file_name,'w')
        fh.write('x,y,z') # write xyz headers 
        [fh.write(','+key) for key in keys] # write attribute headers 
        fh.write('\n') # drop a line 
        for i in range(self.num_elms):
            line = '{:f},{:f},{:f}'.format(x_coords[i],y_coords[i],z_coords[i])
            for key in keys:
                line += ',{:}'.format(self.attr_cache[key][i])
            line += '\n'
            fh.write(line)
        #close file     
        fh.close()
        
    @staticmethod   # find paraview location in windows    
    def findParaview():
        """ Run on windows to find paraview.exe command.
        
        Returns
        -------
        found: bool
            If True the program was able to find where paraview is installed 
            If false then the program could not find paraview in the default 
            install locations. 
        location: str
            if found == True. The string maps to the paraview executable
            if found == False. then 'n/a' is returned.   
        """
        OpSys=platform.system() 
        if OpSys != "Windows":
            raise OSError("This function is only valid on Windows")
        home_dir = os.path.expanduser('~')
        drive_letter = home_dir.split('\\')[0]
        #find paraview in program files?
        path = os.path.join(drive_letter, os.sep,'Program Files')
        contents = os.listdir(path)
        found = False
        for i,pname in enumerate(contents):
            if pname.find("ParaView") != -1:
                para_dir = os.path.join(path,pname)
                found = True
                break
    
        if not found:#try looking in x86 porgram files instead
            path = os.path.join(drive_letter, os.sep,'Program Files (x86)')
            contents = os.listdir(path)
            for i,pname in enumerate(contents):
                if pname.find("ParaView") != -1:
                    para_dir = os.path.join(path,pname)
                    found = True
                    break
        
        if not found:
            return False, 'n/a' 
        else:
            return True, os.path.join(para_dir,'bin','paraview.exe')
        #the string output can be run in the console if it is enclosed in speech
        #marks , ie <"C/program files/ParaView5.X/bin/paraview.exe">
        
    def paraview(self,fname='ResIPy_mesh.vtk',loc=None):
        """ Show mesh in paraview 
        
        Parameters
        -----------
        fname: str,optional
            A vtk file will be written to working directory and then displayed 
            using paraview. 
        loc: str, optional
            Path to paraview excutable, ignored if not using windows. If not provided
            the program will attempt to find where paraview is installed automatically. 
        """
        #formalities 
        if not isinstance(fname,str):
            raise NameError("Excepted string type argument for 'fname'")
        look4 = True
        if loc is not None:
            look4 = False
            if not isinstance(loc,str):
                raise NameError("Excepted string type argument for 'loc'")
                
        self.write_vtk(fname)#write vtk to working directory with all associated attributes
        op_sys = platform.system()#find kernel type
        if op_sys == "Windows":
            if look4: # find where paraview is installed 
                found, cmd_line = self.findParaview()
                cmd_line = '"' + cmd_line + '"' 
                print('paraview location: %s'%cmd_line) # print location to console
                if not found: # raise exception?
                    print("Could not find where paraview is installed")
                    return # exit function 
            else:
                cmd_line = '"' + loc + '"'#use provided location 
            try:
                Popen(cmd_line+' '+fname,shell=True)
            except PermissionError:
                print("Windows has blocked launching paraview, program try running as admin")      
        else:
            Popen(['paraview', fname])
            
    
    def exportTetgenMesh(self,prefix='mesh',zone=None):
        """Export a mesh like the tetgen format for input into E4D. 
        This format is composed of several files. Currently only tested for 
        3D surface array like surveys. 
        
        Parameters
        ----------  
        prefix: string
            Prefix assigned to exported files. Can include path to file. 
        zone: array like, optional
            If using zones in the mesh, this attribute should include the zone
            attribute which is an array of integers identifying the zone 
            associated with each element. By default each element is assigned to 
            zone 1. 
        Notes
        ----------
        Please note routine is experimental and not garanteed to work 100% 
        """
        #error checking / formalities 
        if not isinstance(prefix,str):
            raise NameError('prefix argument is not a string.')
        if zone is not None:
            if len(zone) != self.num_elms:
                raise ValueError('Number of zone array elements given to exportTetgenMesh does not match the number of elements in the mesh.')
        else:
            zone = [1]*self.num_elms # all elements are inside zone 1 
            
        print('## Exporting to tetgen format  ##')# 
            
        print('Computing boundary conditions ...', end ='')# this an intense process:
        #Compute =neighbour and face element matrixes 
        #1) Faces should be counter clockwise 
        #2) Boundary conditions need to be computed, 
        
        #NB: faces on side/bottom of mesh are given a marker of 2. Face on top of mesh denoted 1. 
        #NB: Nodes on side of mesh are given a marker of 2, 1 if on top of the mesh and 0 if inside the mesh 
        con_mat=self.con_matrix
        if self.neigh_matrix is None: # compute neighbour matrix 
            self.computeNeigh() # this will find the element neighbours and the elements which lie on the outside of the mesh! 
            
        neigh = np.array(self.neigh_matrix).T
        out_elem = np.min(neigh, axis=1) == -1 # elements which have a face on the outside of the mesh 
        neigh_trunc = neigh[out_elem]
        con_trunc = np.array(con_mat).T[out_elem]
        
        ## map used in neigh calculation ##         
        map1 = np.array([2,3,4]) -1 
        map2 = np.array([1,4,3]) -1 
        map3 = np.array([1,2,4]) -1 
        map4 = np.array([1,2,3]) -1 
        nmap = np.array([map1,map2,map3,map4])
        node_x = np.array(self.node_x)
        node_y = np.array(self.node_y)
        node_z = np.array(self.node_z)
        
        face_list = []
        ocheck = [] # orientation check, referenced for making face file 
        for i in range(len(neigh_trunc)):
            idx = np.argwhere(neigh_trunc[i]==-1)
            for j in range(len(idx)):
                enodes = con_trunc[i] # element nodes 
                fnodes = con_trunc[i][nmap[idx[j]]][0] # face nodes
                mnode = [enodes[k] not in fnodes for k in range(4)] # find node missing from the face 
                x = np.append(node_x[fnodes],node_x[enodes[mnode]])
                y = np.append(node_y[fnodes],node_y[enodes[mnode]])
                z = np.append(node_z[fnodes],node_z[enodes[mnode]])
                if interp.check_tetra(x,y,z)==2:#points are counter clockwise
                    fnodes_sorted = fnodes
                else: # reorganise so they are counter clockwise
                    fnodes_sorted = np.array((fnodes[0],fnodes[2],fnodes[1]))
                face_list.append(fnodes_sorted)
                
                #now to compute if face boundary condition              
                xm = np.mean(node_x[fnodes]) #determine approx middle of face
                ym = np.mean(node_y[fnodes])
                zm = np.mean(node_z[fnodes])
                dx = abs(max(node_x[fnodes]) - min(node_x[fnodes])) # work out approx face dimensions 
                dy = abs(max(node_y[fnodes]) - min(node_y[fnodes]))
                dz = abs(max(node_z[fnodes]) - min(node_z[fnodes]))
                dd = dx+dy+dz 
                if dy < 1e-16 or dx < 1e-16 : # side of mesh 
                    ocheck.append(0) 
                else: # base or top of mesh 
                    x = np.append(node_x[fnodes_sorted],xm)
                    y = np.append(node_y[fnodes_sorted],ym)
                    z = np.append(node_z[fnodes_sorted],zm+dd)
                    ocheck.append(interp.check_tetra(x,y,z))
                #if ocheck == 0 then the face is on the side of the mesh 
                #if ocheck == 1 then the face is on top of the mesh 
                #if ocheck == 2 then the face is on the bottom of the mesh 
        
        #assign node boundary markers 
        ochecka = np.array(ocheck,dtype=int)
        face_matrix = np.array(face_list)
        #return ochecka, face_matrix
        top_face = ochecka == 1
        other_face = ochecka != 1
        node_bounday_marker = np.zeros(self.num_nodes,dtype=int) # most nodes are inside the mesh 
        top_nodes = np.unique(face_matrix[top_face])
        node_bounday_marker[top_nodes]=1 # these nodes have a value of 1
        
        face_nodes = np.unique(face_matrix[other_face])
        node_bounday_marker[face_nodes]=2 # these nodes have a value of 2 
        print('done')
                
        #output .node file
        print('writing .node file... ',end='')
        fh = open(prefix+'.1.node','w')
        #header line : 
        #<# of points> <dimension (3)> <# of attributes> <boundary markers (0 or 1)>
        fh.write('{:d}\t{:d}\t{:d}\t{:d}\n'.format(self.num_nodes,self.ndims,1,1))
        #all other lines: 
        #<point #> <x> <y> <z> [attributes] [boundary marker]
        for i in range(self.num_nodes):
            line = '{:d}\t{:f}\t{:f}\t{:f}\t{:d}\t{:d}\n'.format((i+1),
                                                        self.node_x[i],
                                                        self.node_y[i],
                                                        self.node_z[i],
                                                        1,
                                                        node_bounday_marker[i])
            fh.write(line)
        fh.write('# exported from meshTools module in ResIPy electrical resistivity processing package')
        fh.close()   
        print('done.')
            
        #output .ele file 
        print('writing .ele file... ',end='')
        fh = open(prefix+'.1.ele','w')
        #First line: <# of tetrahedra> <nodes per tet. (4 or 10)> <region attribute (0 or 1)>
        fh.write('{:d}\t{:d}\t{:d}\n'.format(self.num_elms,self.type2VertsNo(),1))
        #Remaining lines list # of tetrahedra:<tetrahedron #> <node> <node> ... <node> [attribute]
        for i in range(self.num_elms):
            line = '{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\n'.format((i+1),
                                                         self.con_matrix[0][i]+1,#need to add one because of fortran indexing 
                                                         self.con_matrix[1][i]+1,
                                                         self.con_matrix[2][i]+1,
                                                         self.con_matrix[3][i]+1,
                                                         zone[i])
            fh.write(line)
        fh.write('# exported from meshTools module in ResIPy electrical resistivity processing package')
        fh.close()
        print('done.')
        
        #output .trn file 
        print('writing .trn file... ',end='')
        fh = open(prefix+'.trn','w')
        fh.write('0\t0\t0')
        fh.close()
        print('done.')
        
        #write .face file - which describes elements on the outer edges of the mesh
        
        print('Writing .face file... ',end='') 

        fh = open(prefix+'.1.face','w')
        boundary_marker = np.zeros(len(face_list),dtype=int)+2
        boundary_marker[ochecka==1]=1 # set boundary marker to one if top facing element 
        #First line: <# of faces> <boundary marker (0 or 1)>
        fh.write('{:d}\t{:d}\n'.format(len(face_list),1))#header line 
        for i in range(len(face_list)):
            line = '{:d}\t{:d}\t{:d}\t{:d}\t{:d}\n'.format((i+1),
                                                        face_list[i][0]+1,
                                                        face_list[i][1]+1,
                                                        face_list[i][2]+1,
                                                        boundary_marker[i])
            fh.write(line)
            
        
        fh.write('# exported from meshTools module in ResIPy electrical resistivity processing package')    
        fh.close()
        print('done.')
        
        #out .neigh file
        #Here we want look for faces which share with another element 
        #(handled by mesh calc)
            
        neigh_matrix = np.array(self.neigh_matrix).T
        neigh_matrix+=1 # add one for tetgen indexing
        neigh_matrix[neigh_matrix==0]=-1#-1 for face elements 
        
        print('writing .neigh file... ',end='')
        fh = open(prefix+'.1.neigh','w') # write to file 
        fh.write('{:d}\t{:d}\n'.format(self.num_elms,self.type2VertsNo())) # header line         
        elm_id = self.elm_id
        for i in range(self.num_elms):
            line = '{:d}\t{:d}\t{:d}\t{:d}\t{:d}\n'.format(elm_id[i],
                                                            neigh_matrix[i,0],
                                                            neigh_matrix[i,1],
                                                            neigh_matrix[i,2],
                                                            neigh_matrix[i,3])
            fh.write(line)
            
        fh.write('# exported from meshTools module in ResIPy electrical resistivity processing package')    
        fh.close()
        print('done.')
    
#%% import a vtk file 
def vtk_import(file_path='mesh.vtk', order_nodes=True):
    """
    Imports a mesh file into the python workspace, can have triangular, quad or tetraheral shaped elements.
            
    Parameters
    ----------
    file_path : string, optional
        File path to mesh file. Note that a error will occur if the file format is not as expected.
    order_nodes : bool, optional
        Order Nodes if true. Process can be resource intensive though. Reccomended
        if using the mesh for an inversion. 
            
    Returns
    -------
    mesh : class 
        a <ResIPy> mesh class 
    """
    if os.path.getsize(file_path)==0: # So that people dont ask me why you cant read in an empty file, throw up this error. 
        raise ImportError("Provided mesh file is empty! Check that (c)R2/3t code has run correctly!")
    #open the selected file for reading
    fid=open(file_path,'r')
    #print("importing vtk mesh file into python workspace...")
    #read in header info and perform checks to make sure things are as expected
    vtk_ver=fid.readline().strip()#read first line
    if vtk_ver.find('vtk')==-1:
        raise ImportError("Unexpected file type... ")
    elif vtk_ver.find('3.0')==-1:#not the development version for this code
        print("Warning: vtk manipulation code was developed for vtk datafile version 3.0, unexpected behaviour may occur in resulting mesh")
    title=fid.readline().strip()#read line 2
    format_type=fid.readline().strip()#read line 3
    if format_type=='BINARY':
        raise ImportError("expected ASCII type file format, not binary")
    dataset_type=fid.readline().strip().split()#read line 4
    if dataset_type[1]!='UNSTRUCTURED_GRID':
        print("Warning: code is built to parse a vtk 'UNSTRUCTURED_GRID' data type not %s"%dataset_type[1])
    
    #read node data
    #print("importing mesh nodes...")
    node_info=fid.readline().strip().split()#read line 5
    try:
        numnp=int(node_info[1])
    except IndexError:#if we get this then there is a white space between the node info and header lines
        node_info=fid.readline().strip().split()#read line 5
        numnp=int(node_info[1])      
    if numnp == 0: 
        raise ImportError("No nodes in vtk file to import! Aborting... ")
    #now read in node data
    node_x=[0] * numnp #make lists for each of the relevant parameters for each node
    node_y=[0] * numnp
    node_z=[0] * numnp
    node_num=[0] * numnp
    for i in range(numnp):
        try:
            coord_data=fid.readline().strip().split()
            node_x[i] = float(coord_data[0])
            node_y[i] = float(coord_data[1])
            node_z[i] = float(coord_data[2])
        except:# ValueError:
            coord_data=fid.readline()
            node_x[i] = float(coord_data[0:12]) # retrive fixed width columns if cannot parse as split strings
            node_y[i] = float(coord_data[12:24])
            node_z[i] = float(coord_data[24:36])
        node_num[i] = i
    
    #now read in element data
    elm_info=fid.readline().strip().split()#read line with cell data
    try:
        no_elms=int(elm_info[1])
    except IndexError: # quick bug fix
        elm_info=fid.readline().strip().split()#read line with cell data
        no_elms=int(elm_info[1])
    
    if no_elms ==0: 
        raise ImportError("No elements in vtk file to import!")
        
    #read in first element and decide what mesh type it is 
    elm_data=fid.readline().strip().split()
    npere = int(elm_data[0]) # number of vertices per element 
    elm_num = [0]*no_elms
    con_mat = [[0]*no_elms for i in range(npere)]
    for j in range(npere):
        con_mat[j][0] = int(elm_data[j+1])
    
    #read in the rest of the elements 
    for i in range(1,no_elms):
        elm_data=fid.readline().strip().split()
        if int(elm_data[0]) != npere:
            raise ImportError("VTK file contains mixed element types, which are not supported by ResIPy mesh class, aborting...")
        for j in range(npere):
            con_mat[j][i] = int(elm_data[j+1])
        elm_num[i] = i
    
    cell_attr_dump=fid.readlines()#reads the last portion of the file
    #finished reading the file
    
    #find cell types
    for i,line_info in enumerate(cell_attr_dump):
        if line_info.find("CELL_TYPES") == 0:
            #cell_type = [int(k) for k in cell_attr_dump[i+1].strip().split()]
            line_idx = i+1
            cell_type = []
            while len(cell_type) < no_elms:
                cell_type += [float(x) for x in cell_attr_dump[line_idx].split()]
                line_idx += 1
            break
    
    fid.close()
    # read through cell attributes to find the relevant parameter table?
    
    #find scalar values in the vtk file
    num_attr = 0
    attr_dict = {}

    for i,line_info in enumerate(cell_attr_dump):
        if line_info.find("SCALARS") == 0:
            attr_title = line_info.split()[1]
            #check look up table
            if cell_attr_dump[i+1].split()[1] != "default":
                warnings.warn("unrecognised lookup table type")
            line_idx = i+2
            #values=[float(k) for k in cell_attr_dump[i+2].split()]
            array = []
            #while loop to get all scalar values 
            while len(array) < no_elms:
                array += [float(x) for x in cell_attr_dump[line_idx].split()]
                line_idx += 1
            attr_dict[attr_title] = array
            
            if num_attr == 0:# primary attribute defaults to the first attribute found
                parameter_title = attr_title
                values_oi = array   
            num_attr += 1
    
    #put in fail safe if no attributes are found        
    if num_attr == 0:
        values_oi= [0]*no_elms
        parameter_title = "n/a"
        
    #if cell_type[0] == 5 or cell_type[0] == 8 or cell_type[0] == 9: # then its a 2D mesh
    if title == 'Output from cR2' or title == 'Output from R2': # account for the fact the y and z columns should be swapped 
        mesh = Mesh(node_x,#x coordinates of nodes 
                    node_z,#y coordinates of nodes
                    node_y,#z coordinates of nodes 
                    node_num,#node id number 
                    elm_num,#element id number 
                    con_mat,#nodes of element vertices
                    cell_type,#according to vtk format
                    values_oi,#the values of the attributes given to each cell 
                    parameter_title,#what is the attribute? we may use conductivity instead of resistivity for example
                    file_path) #nb: nodes should ordered already  
    else:
        mesh = Mesh(node_x,#x coordinates of nodes 
                    node_y,#y coordinates of nodes
                    node_z,#z coordinates of nodes 
                    node_num,#node id number 
                    elm_num,#element id number 
                    con_mat,#nodes of element vertices
                    cell_type,#according to vtk format
                    values_oi,#the values of the attributes given to each cell 
                    parameter_title,#what is the attribute? we may use conductivity instead of resistivity for example
                    file_path,
                    order_nodes) 
    
    #add attributes / cell parameters 
    for key in attr_dict.keys():
        mesh.attr_cache[key] = attr_dict[key]
    
    #see if sensivity output from R2/R3t is inside the .vtk 
    try:
        if mesh.ndims==2:
            mesh.add_sensitivity(mesh.attr_cache['Sensitivity(log10)'])
        else:
            mesh.add_sensitivity(mesh.attr_cache['Sensitivity_map(log10)'])
    except:
        pass 
    
    mesh.mesh_title = title
    return mesh

#%% import mesh from native .dat format
def dat_import(file_path='mesh.dat', order_nodes=True):
    """ Import R2/cR2/R3t/cR3t .dat kind of mesh. 
    
    Parameters
    ----------
    file_path: str
        Maps to the mesh (.dat) file.
    
    Returns
    -------
    mesh: class
        
    """
    if not isinstance(file_path,str):
        raise TypeError("Expected string type argument for 'file_path'")
    fid=open(file_path,'r')#open file
    #discover if the file is 2 or 3D 
    header = fid.readline().split()
    if len(header)==2: # its a 2D mesh 
        flag_3d = False
        npere = 3
    else:
        flag_3d = True
        npere = int(header[-1])
    numel = int(header[0]) # number of elements 
    numnp = int(header[1]) # number of nodes 
    #now read in element data
    elm_no = [0]*numel # element number / index 
    #allocate nodes
    for i in range(npere):
        exec('node%i = [0]*%i'%(i,numel))
    node_map =[[0]*numel for i in range(npere)] 
    #np.array([[0]*numel]*npere,dtype=int)
    zone = [0]*numel # mesh zone 
    for i in range(numel):
        line = fid.readline().split()#read in line data
        elm_no[i] = int(line[0])
        zone[i] = int(line[-1])
        ref=1 # read through each node index 
        if i==0 and not flag_3d and len(line)==7: # we have a special case of a quad mesh .dat file
            npere=4
        for j in range(npere):
            #exec('node%i[i] = int(line[ref])-1'%j)
            node_map[j][i]=int(line[ref])-1
            ref += 1

    #read in nodes 
    node_x = [0]*numnp
    node_y = [0]*numnp
    node_z = [0]*numnp
    node_id = [0]*numnp
    for i in range(numnp):
        line = fid.readline().split()#read in line data
        while len(line)==0: # then we have a space between elements and nodes 
            line = fid.readline().split()#read in line data
        node_id[i] = int(line[0])
        node_x[i] = float(line[1])
        if flag_3d:
            node_y[i] = float(line[2])
            node_z[i] = float(line[3])
        else:
            node_y[i] = float(0)
            node_z[i] = float(line[2])

    
    #probe vtk cell type
    if flag_3d:
        if npere == 4:
            cell_type = 10
        elif npere == 6:
            cell_type = 13
    else:
        if npere == 4:
            cell_type = 9
        elif npere == 3:
            cell_type = 5
            
    #iniate mesh class 
    mesh = Mesh(node_x = node_x,#x coordinates of nodes 
                node_y = node_y,#y coordinates of nodes
                node_z = node_z,#z coordinates of nodes 
                node_id= node_id,#node id number 
                elm_id=elm_no,#element id number 
                node_data=node_map,#nodes of element vertices
                cell_type = [cell_type],#according to vtk format
                cell_attributes = zone,#the values of the attributes given to each cell, we dont have any yet 
                atribute_title='zone', #what is the attribute? we may use conductivity instead of resistivity for example
                order_nodes = order_nodes)
    
    mesh.add_attribute(zone,'zone')
    
    return mesh 
        
           
#%% Read in E4D / tetgen mesh
def tetgen_import(file_path, order_nodes=True):
    """Import Tetgen mesh into ResIPy. This isa little different from other 
    imports as the mesh is described by several files. From meshTools' perspective
    What is needed is the node(.node) file and element (.ele) files, which 
    describe the coordinates of mesh nodes and the connection matrix. 
    
    Parameters
    ----------
    file_path: str
        Maps to the mesh .node file. The program will automatically find the 
        corresponding .ele file in the same directory. 
    
    Returns
    -------
    mesh: class    
    """
    fh = open(file_path,'r') # open file for read in 
    header = fh.readline() # header line
    numnp = int(header.split()[0]) # number of nodes 
    node_x = [0]*numnp # preallocate array likes for node coordinates 
    node_y = [0]*numnp
    node_z = [0]*numnp
    node_id = [0]*numnp
    for i in range(numnp): # read through each line and get node coordinates 
        line = fh.readline().split()
        node_id[i] = int(line[0])
        node_x[i] = float(line[1])
        node_y[i] = float(line[2])
        node_z[i] = float(line[3])
    fh.close() #close file 
    
    #try and find the .trn file? - this is an output from E4D which describes 
    # a mesh translation applied to improve computational accuracy 
    file_path_trn =[file_path.replace('.node','.trn'),
                    file_path.replace('.1.node','.trn')]# 2 likely options for the way the file is named 
    for i in range(2):
        if os.path.exists(file_path_trn[i]):
            fh = open(file_path_trn[i],'r')
            line = fh.readline().split()
            fh.close()
            delta_x = float(line[0])
            delta_y = float(line[1])
            delta_z = float(line[2])
            break # break loop if found the file 
        else:
            delta_x = 0
            delta_y = 0
            delta_z = 0
    
    node_x = np.array(node_x) + delta_x
    node_y = np.array(node_y) + delta_y
    node_z = np.array(node_z) + delta_z
    
    #next we need the connection matrix to describe each of the tetrahedral e
    #elements
    file_path2 = file_path.replace('.node','.ele')
    fh = open(file_path2,'r')# read in element file  
    header = fh.readline() # read in header line 
    numel = int(header.split()[0]) # number of elements 
    #meshes so its always going to be 4. 
    
    
    node_map = ([0]*numel,[0]*numel,[0]*numel,[0]*numel)
    #np.array([[0]*numel]*npere,dtype=int) # connection matrix mapping elements onto nodes 
    elm_no = [0]*numel # element number / index 
    zone = [0]*numel # mesh zone 
    
    for i in range(numel):
        line = fh.readline().split()#read in line data
        elm_no[i] = int(line[0])
        zone[i] = int(line[-1])
        node_map[0][i]=int(line[1])-1
        node_map[1][i]=int(line[2])-1
        node_map[2][i]=int(line[3])-1
        node_map[3][i]=int(line[4])-1
    
    #create mesh instance 
    mesh = Mesh(node_x = node_x,#x coordinates of nodes 
                node_y = node_y,#y coordinates of nodes
                node_z = node_z,#z coordinates of nodes 
                node_id= node_id,#node id number 
                elm_id=elm_no,#element id number 
                node_data=node_map,#nodes of element vertices
                cell_type = [10],#according to vtk format
                cell_attributes = zone,#the values of the attributes given to each cell, we dont have any yet 
                original_file_path = file_path,
                atribute_title='zone',#what is the attribute? 
                order_nodes = order_nodes)
    
    mesh.add_attribute(zone,'zone')
    
    return mesh
        
        
#%% build a quad mesh        
def quad_mesh(elec_x, elec_z, elec_type = None, elemx=4, xgf=1.5, zf=1.1, zgf=1.25, fmd=None, pad=2, 
              surface_x=None,surface_z=None):
    """Creates a quaderlateral mesh given the electrode x and y positions.
            
    Parameters
    ----------
    elec_x : list, np array
        Electrode x coordinates 
    elec_z : list, np array
        Electrode y coordinates
    elec_type: list, optional
        strings, where 'electrode' is a surface electrode; 'buried' is a buried electrode
    elemy : int, optional
        Number of elements in the fine y region
    zf : float, optional
         Z factor multiplier in the fine zone.
    zgf : float, optional
         Z factor multiplier in the coarse zone.
    fmd : float (m), optional 
         Fine mesh region depth specifies as positive number (if None, half survey width is used).
    pad : int, optional
         X padding outside the fine area (tipicaly twice the number of elements between electrodes).
    surface_x: array like, optional
        Default is None. x coordinates of extra surface topography points for the generation of topography in the quad mesh
    surface_z: array like, optional
        Default is None. z coordinates of extra surface topography points for the generation of topography in the quad mesh. Note
        an error will be returned if len(surface_x) != len(surface_z)
        
            
    Returns
    -------
    Mesh : class
        Mesh object 
    meshx : numpy.array
        Mesh x locations for R2in file.
    meshz : numpy.array
        Mesh  locations for R2in file (ie node depths).
    topo : numpy.array
        Topography for R2in file.
    elec_node : numpy.array
        x columns where the electrodes are. 
    """        
    if surface_x is None or surface_z is None:
        surface_x = np.array([])
        surface_z = np.array([])
    else: #if surface_x != None and surface_y != None:
        if len(surface_x) != len(surface_z):
            raise Exception("The length of the surface_x argument does not match the surface_y argument, both need to be arrays of the same length.")
        surface_x = np.array(surface_x)
        surface_z = np.array(surface_z)
    
    bh_flag = False
    #determine the relevant node ordering for the surface electrodes? 
    if elec_type is not None:
        if not isinstance(elec_type,list):
            raise TypeError("electrode_type variable is given but expected a list type argument")
        if len(elec_type) != len(elec_x):
            raise ValueError("mis-match in length between the electrode type and number of electrodes")
        
        surface_idx=[]#surface electrode index
        bur_idx=[]#buried electrode index 
        for i,key in enumerate(elec_type):
            if key == 'electrode': surface_idx.append(i)
            if key == 'buried': 
                bur_idx.append(i)
                bh_flag = True
        
        if len(surface_idx)>0:# then surface electrodes are present
            Ex = np.array(elec_x)[surface_idx]
            Ey = np.array(elec_z)[surface_idx]
        elif len(surface_idx) == 0 and len(surface_x) > 0:
            #case where you have surface topography but no surface electrodes 
            Ex = np.array(surface_x)
            Ey = np.array(surface_z)
            elec = np.c_[Ex, Ey]
        elif len(surface_idx)== 0:
            #fail safe if no surface electrodes are present to generate surface topography 
            Ex = np.array([elec_x[np.argmin(elec_x)], elec_x[np.argmax(elec_x)]])
            Ey = np.array([elec_z[np.argmax(elec_z)], elec_z[np.argmax(elec_z)]])
    else:
        pass
        #elec = np.c_[elec_x,elec_y]
    elec = np.c_[elec_x,elec_z]    
    
    if bh_flag:
        bh = np.c_[np.array(elec_x)[bur_idx],np.array(elec_z)[bur_idx]]
        
    pad = pad # number of padding on both side (as a multiplier of the nb of nodes between electrodes)
    # create meshx
    meshx = np.array([])
    elecXsorted=np.sort(elec[:,0]) # sort x coordinates of the electrodes 
    for i in range(len(elec)-1):
        elec1 = elecXsorted[i]#elec[i,0]
        elec2 = elecXsorted[i+1]#elec[i+1,0]
        espacing = np.abs(elec1-elec2)
#        dx = espacing/elemx # we ask for elemx nodes between electrodes
        if espacing == 0: # then we probably are probing 2 borehole electrodes
            espacing = np.abs(elec[i,1] - elec[i+1,1])
        if i == 0:
            xx2 = np.linspace(elec1-espacing, elec1, elemx, endpoint=False)
            xx3 = np.ones(elemx*pad)*elec1-espacing
            dxx = espacing
            for j in range(1,elemx*pad): # padding
                xx3[j] = xx3[j-1]-dxx*xgf
                dxx = dxx*xgf
            meshx = np.r_[meshx, xx3[::-1], xx2[1:]]
        xx = np.linspace(elec1, elec2, elemx, endpoint=False)
        meshx = np.r_[meshx, xx]
        if i == len(elec)-2:
            xx2 = np.linspace(elec2, elec2+espacing, elemx, endpoint=False)
            xx3 = np.ones(elemx*pad)*elec2+espacing
            dxx = espacing
            for j in range(1,elemx*pad):
                xx3[j] = xx3[j-1]+dxx*xgf
                dxx = dxx*xgf
            meshx = np.r_[meshx, xx2, xx3]
            
    # create meshz
    if fmd is None:
        if bh_flag:
            fmd = abs(min(elec_z))
        else:
            fmd = (max(elec[:,0]) - min(elec[:,0]))/2
        
#    dyy = espacing/(elemx*4)
    meshz = [0]
    dyy = 0.05
    for i in range(100):
        meshz.append(meshz[-1]+dyy*zf)
        dyy = dyy*zf
        if meshz[-1] > fmd:
            break
    elemy = len(meshz)
    elemy2 = int(elemy/2)
    yy = np.ones(elemy2)*meshz[-1]
    for i in range(1, elemy2):
        yy[i] = yy[i-1]+dyy*zgf
        dyy = dyy*zgf
    meshz = np.r_[meshz, yy[1:]]
    
    #insert borehole electrodes? if we have boreholes / buried electrodes 
    if bh_flag:
        meshx = np.unique(np.append(meshx,bh[:,0]))
        
    # create topo
    if bh_flag: # only use surface electrodes to make the topography if buried electrodes present
        X = np.append(Ex, surface_x) 
        Y = np.append(Ey, surface_z)
        idx = np.argsort(X)
        topo = np.interp(meshx, X[idx], Y[idx])
    else: # all electrodes are assumed to be on the surface 
        X = np.append(elec[:,0], surface_x)
        Y = np.append(elec[:,1], surface_z)
        idx = np.argsort(X)
        topo = np.interp(meshx, X[idx], Y[idx])
    
    if bh_flag:
        #insert y values of boreholes, normalised to topography
        norm_bhy = np.interp(bh[:,0], elec[:,0], elec[:,1]) - bh[:,1]
        meshz = np.unique(np.append(meshz,norm_bhy))
    
    ###
    #find the columns relating to the electrode nodes? 
    temp_x = meshx.tolist()
    temp_y = meshz.tolist()
    elec_node_x=[temp_x.index(elec_x[i])+1 for i in range(len(elec_x))]#add 1 because of indexing in R2. 
    
    elec_node_y = np.ones((len(elec_z),1),dtype=int)#by default electrodes are at the surface
    #this ugly looking if - for loop thing finds the y values for borehole meshes. 
    if bh_flag:
        count = 0 
        try:            
            for i,key in enumerate(elec_type):
                if key == 'buried':
                    elec_node_y[i] = temp_y.index(norm_bhy[count]) + 1
                    #print(elec_node_y[i])
                    count += 1
        except ValueError:
            raise Exception("There was a problem indexing the meshz values for electrode number %i"%i)
                
    elec_node = [elec_node_x,elec_node_y.T.tolist()[0]]
    
    #print some warnings for debugging 
    if len(topo)!=len(meshx):
        warnings.warn("Topography vector and x coordinate arrays not the same length! ")
    elif len(elec_node_x)!=len(elec_x):
        warnings.warn("Electrode node vector and number of electrodes mismatch! ")
     
    # what is the number of regions? (elements)
    no_elms=(len(meshx)-1)*(len(meshz)-1)
    no_nodes=len(meshx)*len(meshz)
    
    # compute node mappins (connection matrix)
    y_dim=len(meshz)
    fnl_node=no_nodes-1
    
    node_mappins=(np.arange(0,fnl_node-y_dim), # top left
                  np.arange(1,fnl_node-y_dim+1), # bottom left
                  np.arange(y_dim+1,fnl_node+1), # bottom right
                  np.arange(y_dim,fnl_node)) # top right
    
    del_idx = np.arange(y_dim-1,len(node_mappins[0]),y_dim)#the above has too many indexes at the changeover of columns so some need deleting
    
    node_mappins = [list(np.delete(node_mappins[i],del_idx)) for i in range(4)]#delete excess node placements
    #compute node x and y  (and z)
    node_x,node_z = np.meshgrid(meshx,meshz)
    #account for topography in the y direction 
    node_z = [topo-node_z[i,:] for i in range(y_dim)]#list comprehension to add topography to the mesh, (could use numpy here??)
    node_z = np.array(node_z).flatten(order='F')
    node_x = node_x.flatten(order='F')
    node_y = np.array([0]*len(node_x))
       
    #make mesh class    
    mesh = Mesh(node_x,
                node_y,
                node_z,
                list(np.arange(0,no_nodes)),
                list(np.arange(0,no_elms)),
                node_mappins,
                [9],
                [1]*no_elms,
                'zones',
                order_nodes=False)
    
    #find the node which the electrodes are actually on in terms of the mesh. 
    node_in_mesh = [0]*len(elec_x)
    for i in range(len(elec_x)):
        sq_dist = (node_x - elec_x[i])**2 + (node_z - elec_z[i])**2 # find the minimum square distance
        
        node_in_mesh[i] = np.argmin(sq_dist) # min distance should be zero, ie. the node index.
    mesh.add_e_nodes(node_in_mesh) # add nodes to the mesh class
    
    # point at the surface
    xsurf = []
    zsurf = []
    for x, z, t in zip(elec_x, elec_z, elec_type):
        if t == 'electrode': # surface electrode
            xsurf.append(x)
            zsurf.append(z)
    if surface_x is not None:
        xsurf = xsurf + list(surface_x)
        zsurf = zsurf + list(surface_z)
    surfacePoints = np.array([xsurf, zsurf]).T
    isort = np.argsort(xsurf)
    mesh.surface = surfacePoints[isort, :]

    # OVERWRITE elec_node using node_in_mesh that doesn't assume surface
    # electrode but snap them to node based on distance
    elec_node = np.c_[np.arange(len(node_in_mesh))+1, node_in_mesh]

    return mesh, meshx, meshz, topo, elec_node

#%% build a triangle mesh - using the gmsh wrapper
def tri_mesh(elec_x, elec_z, elec_type=None, geom_input=None, keep_files=True, 
             show_output=True, path='exe', dump=print, whole_space=False, 
             **kwargs):
    """ Generates a triangular mesh for r2. Returns mesh class ...
    this function expects the current working directory has path: exe/gmsh.exe.
    Uses gmsh version 3.0.6.
            
    Parameters
    ---------- 
    elec_x: array like
        electrode x coordinates 
    elec_z: array like 
        electrode y coordinates 
    elec_type: list of strings, optional
        List should be the same length as the electrode coordinate argument. Each entry in
        the list references the type of electrode: 
        - 'electrode' = surface electrode coordinate, will be used to construct the topography in the mesh
        - 'buried' = buried electrode, placed the mesh surface
        - 'borehole' = borehole electrode, electrodes will be placed in the mesh with a line connecting them. 
        borehole numbering starts at 1 and ascends numerically by 1.  
    geom_input : dict, optional
        Allows for further customisation of the 2D mesh, its a
        dictionary contianing surface topography, polygons and boundaries 
    keep_files : boolean, optional
        `True` if the gmsh input and output file is to be stored in the exe directory.
    show_ouput : boolean, optional
        `True` if gmsh output is to be printed to console. 
    path : string, optional
        Path to exe folder (leave default unless you know what you are doing).
    whole_space: boolean, optional
        flag for if the problem should be treated as a whole space porblem, in which case 
        electrode type is ingored and all electrodes are buried in the middle of a large mesh. 
    dump : function, optional
        Function to which pass the output during mesh generation. `print()` is
        the default.
    **kwargs : optional
        Key word arguments to be passed to genGeoFile. 
            
    Returns
    -------
    mesh: class
        <ResIPy> mesh class
        
    Notes
    -----
    geom_input format:
        the code will cycle through numerically ordered keys (strings referencing objects in a dictionary"),
        currently the code expects a 'surface' and 'electrode' key for surface points and electrodes.
        the first borehole string should be given the key 'borehole1' and so on. The code stops
        searching for more keys when it cant find the next numeric key. Same concept goes for adding boundaries
        and polygons to the mesh. See below example:
            
            geom_input = {'surface': [surf_x,surf_z],
              'boundary1':[bound1x,bound1y],
              'polygon1':[poly1x,poly1y]} 
            
    electrodes and electrode_type (if not None) format: 
        
            electrodes = [[x1,x2,x3,...],[y1,y2,y3,...]]
            electrode_type = ['electrode','electrode','buried',...]
        
        like with geom_input, boreholes should be labelled borehole1, borehole2 and so on.
        The code will cycle through each borehole and internally sort them and add them to 
        the mesh. 
        
    The code expects that all polygons, boundaries and electrodes fall within x values 
    of the actaul survey area. So make sure your topography / surface electrode points cover 
    the area you are surveying, otherwise some funky errors will occur in the mesh. 
    """
    #check directories 
    if path == "exe":
        ewd = os.path.join(
                os.path.dirname(os.path.realpath(__file__)),
                path)
        print('exe dir = ' + ewd) #ewd - exe working directory 
    else:
        ewd = path
        # else its assumed a custom directory has been given to the gmsh.exe
    
    if not os.path.isfile(os.path.join(ewd,'gmsh.exe')) and not os.path.isfile(os.path.join(ewd,'gmsh_linux')):
        raise Exception("No gmsh executable exists in the exe directory!")
    
    #make .geo file
    file_name="mesh"
    if not whole_space:#by default create survey with topography 
        node_pos = gw.genGeoFile([elec_x,elec_z], elec_type, geom_input,
                             file_path=file_name,**kwargs)
    elif whole_space:
        print("Whole space problem")
        node_pos = gw.gen_2d_whole_space([elec_x,elec_z], geom_input = geom_input, 
                                         file_path=file_name,**kwargs)    
    
    # handling gmsh
    if platform.system() == "Windows":#command line input will vary slighty by system 
#        cmd_line = ewd+'\gmsh.exe '+file_name+'.geo -2'
        cmd_line = os.path.join(ewd,'gmsh.exe')+' '+file_name+'.geo -2'
    elif platform.system() == 'Darwin':
            winePath = []
            wine_path = Popen(['which', 'wine'], stdout=PIPE, shell=False, universal_newlines=True)#.communicate()[0]
            for stdout_line in iter(wine_path.stdout.readline, ''):
                winePath.append(stdout_line)
            if winePath != []:
                cmd_line = ['%s' % (winePath[0].strip('\n')), ewd+'/gmsh.exe', file_name+'.geo', '-2']
            else:
                cmd_line = ['/usr/local/bin/wine', ewd+'/gmsh.exe', file_name+'.geo', '-2']
    else:
        if os.path.isfile(os.path.join(ewd,'gmsh_linux')):
            cmd_line = [ewd + '/gmsh_linux', file_name + '.geo', '-2'] # using linux version if avialable (can be more performant)
        else: # fall back to wine
            cmd_line = ['wine',ewd+'/gmsh.exe', file_name+'.geo', '-2']

    if show_output: 
        p = Popen(cmd_line, stdout=PIPE, stderr=PIPE, shell=False)#run gmsh with ouput displayed in console
        while p.poll() is None:
            line = p.stdout.readline().rstrip()
            if line.decode('utf-8') != '':
                dump(line.decode('utf-8'))
    else:
        p = Popen(cmd_line, stdout=PIPE, stderr=PIPE, shell=False)
        p.communicate() # wait to finish
        
    #convert into mesh.dat 
    mesh_info = gw.msh_parse(file_path = file_name+'.msh') # read in mesh file
    
    # merge fine with coarse regions
    regions = np.array(mesh_info['parameters'])
    for reg in np.unique(regions)[1:]:
        ie = regions == reg
        regions[ie] = reg - 1
    
    mesh = Mesh(mesh_info['node_x'], # convert output of parser into an object
                mesh_info['node_y'],
                mesh_info['node_z'],
                mesh_info['node_id'],
                mesh_info['elm_id'],
                mesh_info['node_data'],
                mesh_info['cell_type'],
                regions, # this is the region
                mesh_info['parameter_title'], # and it's label
                mesh_info['original_file_path'])
    
    if keep_files is False: 
        os.remove(file_name+".geo")
        os.remove(file_name+".msh")

    mesh.add_e_nodes(node_pos-1)#in python indexing starts at 0, in gmsh it starts at 1 
    
    # point at the surface
    xsurf = []
    zsurf = []
    # add remote if any
    if elec_type is not None:
        iremote = np.array([a == 'remote' for a in elec_type])
        mesh.iremote = iremote
    
        for x, z, t in zip(elec_x, elec_z, elec_type):
            if t == 'electrode': # surface electrode
                xsurf.append(x)
                zsurf.append(z)
    
    if geom_input is not None:             
        if 'surface' in geom_input.keys():
            xsurf = xsurf + list(geom_input['surface'][0])
            zsurf = zsurf + list(geom_input['surface'][1])
        surfacePoints = np.array([xsurf, zsurf]).T
        isort = np.argsort(xsurf)
        mesh.surface = surfacePoints[isort, :]

    return mesh


#%% 3D tetrahedral mesh 
def tetra_mesh(elec_x,elec_y,elec_z=None, elec_type = None, keep_files=True, interp_method = 'triangulate',
               surface_refinement = None, mesh_refinement = None,show_output=True, 
               path='exe', dump=print, whole_space=False, padding=20, ncores=2,
               search_radius = 10, **kwargs):
    """ Generates a tetrahedral mesh for R3t (with topography). returns mesh3d.dat 
    in the working directory. This function expects the current working directory 
    has path: exe/gmsh.exe.
    
    Uses post processing after mesh generation to super impose topography on to 
    a flat 3D tetrahedral mesh. 
            
    Parameters
    ---------- 
    elec_x: array like
        electrode x coordinates 
    elec_y: array like 
        electrode y coordinates 
    elec_z: array like 
        electrode z coordinates 
    elec_type: list of strings, optional
        Defines if electrodes are buried or not.   
    keep_files : boolean, optional
        `True` if the gmsh input and output file is to be stored in the working directory.
    interp_method: string, default ='bilinear' optional
        Interpolation method to translate mesh nodes in the z direction. In other words the method in which topography 
        is appended to the mesh. Here the topography is added to the mesh in post processing. 
        The options documented in the notes. 
        if == 'idw': then provide search_radius.  
    surface_refinement : np.array, optional 
        Numpy array of shape (3,n), should follow the format np.array([x1,x2,x3,...],[y1,y2,y3,...],[z1,z2,z3,...]).
        Allows for extra refinement for the top surface of the mesh. The points are not added to the mesh, but 
        considered in post processing in order to super impose topography on the mesh. 
    mesh_refinement : np.array, optional
        Coming soon ... Not yet implimented.
    show_ouput : boolean, optional
        `True` if gmsh output is to be printed to console. 
    path : string, optional
        Path to exe folder (leave default unless you know what you are doing).
    whole_space: boolean, optional
        flag for if the problem should be treated as a whole space porblem, in which case 
        electrode type is ingored and all electrodes are buried in the middle of a large mesh. 
    dump : function, optional
        Function to which pass the output during mesh generation. `print()` is
        the default.
    padding : float, optional
        amount of % the fine mesh region will be extended beyond the extent of the electrode positions
    search_radius: float, None, optional
        Defines search radius used in the inverse distance weighting interpolation. 
        If None then no search radius will be used and all points will be considered in the interpolation. 
    **kwargs : optional
        Key word arguments to be passed to box_3d. 
            
    Returns
    ---------- 
    mesh3d: class
 
    Notes 
    ---------- 
    Possible arguments for interp_method: 
        'bilinear' : 4 known points are used to compute the equation of a plane 
                    in which the interpolated point lies. This method is reccomended 
                    if elevation data is organised in a regular grid. 
        'idw' : Inverse distance weighting, interpolated points are assigned a 
                    a z value which is a weighted average of all points in the 
                    search raduis. This method works best for gentle topography. 
        'nearest' : Nearest neighbour interpolation. Z value at the interpolated 
                    point takes on the same Z value as the closest known point. 
                    This method can work well for dense elevation data, say in 
                    the case using a point cloud from a digital elevation model. 
                    This method is the most computationally cheap to run. 
        'spline' : Like bilinear interpolation, a function is fitted to a 
                    quadlerateral of known points encompassing the interpolated
                    Z value. The computed function is a higher order than linear 
                    interpolation though. Works best for irregular spaced 
                    elevation data. 
        'triangulate': Triangulation method, best for complicated topographies
                    where an irregular grid of known topography points is not 
                    avialable. This is the default. 
        None : No interpolation method is used to interpolated topography on to 
                    mesh, hence a flat mesh is returned. 
    """
    #formalities 
    if len(elec_x) != len(elec_y):
        raise ValueError("mismatch in electrode x and y vector length, they should be equal")
    elif elec_type is not None:
        if len(elec_x) != len(elec_z):
            raise ValueError("mismatch in electrode x and z vector length, they should be equal")
        check = 0
        for i, key in enumerate(elec_type):
            if key=='surface': check+=1
        if len(elec_type)==check:
            print("all electrodes are surface electrodes, ignoring the electrode type")
            elec_type = None
    avail_methods = ['bilinear','idw','nearest','spline','triangulate',None]
    if interp_method not in avail_methods:
        raise NameError("'%s' is an unrecognised interpretation method"%interp_method)
            
    if surface_refinement is not None:
        surf_x = surface_refinement[:,0]
        surf_y = surface_refinement[:,1]
        surf_z = surface_refinement[:,2]
    else:
        surf_x = []
        surf_y = []
        surf_z = []
    
    rem_elec_idx = []
    if elec_type is not None:
        if not isinstance(elec_type,list):
            raise TypeError("'elec_type' argument should be of type 'list', got type %s"%str(type(elec_type)))
        elif len(elec_type) != len(elec_x):
            raise ValueError("mismatch in elec_type vector and elec_x vector lengths")
        surf_elec_x = []
        surf_elec_y = []
        surf_elec_z = []
        bur_elec_x = []
        bur_elec_y = []
        bur_elec_z = []
        surf_elec_idx = []
        bur_elec_idx = []
        
        if elec_z is None:
            elec_z = np.zeros_like(elec_x)
        for i, key in enumerate(elec_type):
            if key == 'buried':
                bur_elec_x.append(elec_x[i])
                bur_elec_y.append(elec_y[i])
                bur_elec_z.append(elec_z[i])
                bur_elec_idx.append(i)
            if key == 'surface' or key=='electrode':
                surf_elec_x.append(elec_x[i])
                surf_elec_y.append(elec_y[i])
                surf_elec_z.append(elec_z[i])
                surf_elec_idx.append(i)
            if key == 'remote':
                rem_elec_idx.append(i)
        #interpolate in order to normalise buried electrode elevations to 0
        x_interp = np.append(surf_elec_x,surf_x)#parameters to be interpolated with
        y_interp = np.append(surf_elec_y,surf_y)
        z_interp = np.append(surf_elec_z,surf_z)
        
        if len(bur_elec_x)>0: #if we have buried electrodes normalise their elevation to as if they are on a flat surface
            print('found buried electrodes')
            if interp_method == 'idw': 
                bur_elec_z = interp.idw(bur_elec_x, bur_elec_y, x_interp, y_interp, z_interp,radius=search_radius)# use inverse distance weighting
            elif interp_method == 'bilinear' or interp_method == None: # still need to normalise electrode depths if we want a flat mesh, so use biliner interpolation instead
                bur_elec_z = interp.interp2d(bur_elec_x, bur_elec_y, x_interp, y_interp, z_interp)
            elif interp_method == 'nearest':
                bur_elec_z = interp.nearest(bur_elec_x, bur_elec_y, x_interp, y_interp, z_interp)
            elif interp_method == 'spline':
                bur_elec_z = interp.interp2d(bur_elec_x, bur_elec_y, x_interp, y_interp, z_interp,method='spline')
            elif interp_method == 'triangulate':
                bur_elec_z = interp.triangulate(bur_elec_x, bur_elec_y, x_interp, y_interp, z_interp)
            elif interp_method is None:
                bur_elec_z = np.zeros_like(bur_elec_idx)
                
        elec_z = np.array(elec_z) 
        elec_z[surf_elec_idx] = 0
        elec_z[bur_elec_idx] = elec_z[bur_elec_idx] - bur_elec_z # normalise to zero surface 
        
    else:
        surf_elec_x = elec_x 
        surf_elec_y = elec_y 
        if elec_z is None:
            surf_elec_z = np.zeros_like(elec_x)
            elec_z = np.zeros_like(elec_x)
        else:
            surf_elec_z = elec_z.copy()
            elec_z = np.array(elec_z) - np.array(elec_z)#normalise elec_z
            
    #check if remeote electrodes present, and remove them 
    if len(rem_elec_idx)>0:
        elec_x = np.delete(elec_x,rem_elec_idx)
        elec_y = np.delete(elec_y,rem_elec_idx)
        elec_z = np.delete(elec_z,rem_elec_idx)
         
    #check directories 
    if path == "exe":
        ewd = os.path.join(
                os.path.dirname(os.path.realpath(__file__)),
                path)
    else:
        ewd = path # points to the location of the .exe 
        # else its assumed a custom directory has been given to the gmsh.exe 
    
    if not os.path.isfile(os.path.join(ewd,'gmsh.exe')) and not os.path.isfile(os.path.join(ewd,'gmsh_linux')):
        raise Exception("No gmsh executable exists in the exe directory!")
    
    #make .geo file
    file_name="mesh3d"
    if whole_space:#by default create survey with topography 
        print("Whole space problem")
        raise Exception("Sorry whole space 3D problems are not implimented yet")
        
    else:
        node_pos = gw.box_3d([elec_x,elec_y,elec_z], file_path=file_name, **kwargs)
            
    # handling gmsh
    if platform.system() == "Windows":#command line input will vary slighty by system 
#        cmd_line = ewd+'\gmsh.exe '+file_name+'.geo -3'
        cmd_line = os.path.join(ewd,'gmsh.exe')+' '+file_name+'.geo -3 -nt %i'%ncores 
    elif platform.system() == 'Darwin':
            winePath = []
            wine_path = Popen(['which', 'wine'], stdout=PIPE, shell=False, universal_newlines=True)#.communicate()[0]
            for stdout_line in iter(wine_path.stdout.readline, ''):
                winePath.append(stdout_line)
            if winePath != []:
                cmd_line = ['%s' % (winePath[0].strip('\n')), ewd+'/gmsh.exe', file_name+'.geo', '-3', 'nt','%i'%ncores]
            else:
                cmd_line = ['/usr/local/bin/wine', ewd+'/gmsh.exe', file_name+'.geo', '-3', 'nt','%i'%ncores]
    else:
        if os.path.isfile(os.path.join(ewd,'gmsh_linux')): # if linux gmsh is present
            cmd_line = [ewd+'/gmsh_linux', file_name+'.geo', '-3', 'nt','%i'%ncores]
        else: # fallback on wine
            cmd_line = ['wine',ewd+'/gmsh.exe', file_name+'.geo', '-3', 'nt','%i'%ncores]
        
    if show_output: 
        # try:
            # p = Popen(cmd_line, stdout=PIPE, shell=False)#run gmsh with ouput displayed in console
        # except: # hotfix to deal with failing commits on gitlab's server. 
            # cmd_line = ['wine', ewd + '/gmsh.exe', file_name + '.geo', '-3', 'nt', '%i'%ncores] # use .exe through wine instead
        p = Popen(cmd_line, stdout=PIPE, shell=False)
        while p.poll() is None:
            line = p.stdout.readline().rstrip()
            dump(line.decode('utf-8'))
    else:
        p = Popen(cmd_line, stdout=PIPE, stderr=PIPE, shell=False)
        p.communicate() # wait to finish

        
    #convert into mesh.dat
    mesh_info = gw.msh_parse(file_path = file_name+'.msh') # read in 3D mesh file
    
    # merge fine with coarse regions
    regions = np.array(mesh_info['parameters'])
    for reg in np.unique(regions)[1:]:
        ie = regions == reg
        regions[ie] = reg - 1
    
    mesh = Mesh(mesh_info['node_x'], # convert output of parser into an object
                mesh_info['node_y'],
                mesh_info['node_z'],
                mesh_info['node_id'],
                mesh_info['elm_id'],
                mesh_info['node_data'],
                mesh_info['cell_type'],
                regions, # this is used as cell_attributes
                mesh_info['parameter_title'],
                mesh_info['original_file_path'])

    #mesh.write_dat(file_path='mesh.dat') # write mesh.dat - disabled as handled higher up in the R2 class 
    node_x = np.array(mesh.node_x)
    node_y = np.array(mesh.node_y)
    
    if keep_files is False: 
        os.remove(file_name+".geo");os.remove(file_name+".msh")
        
    print('interpolating topography onto mesh using %s interpolation...'%interp_method, end='')
    
    x_interp = np.append(surf_elec_x,surf_x)#parameters to be interpolated with
    y_interp = np.append(surf_elec_y,surf_y)
    z_interp = np.append(surf_elec_z,surf_z)

    #using home grown functions to interpolate / extrapolate topography on mesh
    if interp_method == 'idw': 
        nodez = interp.idw(node_x, node_y, x_interp, y_interp, z_interp,radius=search_radius)# use inverse distance weighting
    elif interp_method == 'bilinear':# interpolate on a irregular grid, extrapolates the unknown coordinates
        nodez = interp.interp2d(node_x, node_y, x_interp, y_interp, z_interp)
    elif interp_method == 'nearest':
        nodez = interp.nearest(node_x, node_y, x_interp, y_interp, z_interp)
    elif interp_method == 'spline':
        nodez = interp.interp2d(node_x, node_y, x_interp, y_interp, z_interp,method='spline')
    elif interp_method == 'triangulate':
        nodez = interp.triangulate(node_x, node_y, x_interp, y_interp, z_interp)
    elif interp_method == None:
        nodez = np.zeros_like(node_x,dtype=float)
    print('done')
    
    mesh.node_z = np.array(mesh.node_z) + nodez
    #need to recompute cell centres as well as they will have changed. 
    mesh.cellCentres()

    #check if remeote electrodes present, and insert them into the node position array
    if len(rem_elec_idx)>0: 
        rem_node_bool = min(node_x) == node_x #& (min(node_y) == node_y) & (min(node_z) == node_z)
        rem_node = np.argwhere(rem_node_bool == True)
        print(rem_elec_idx)
        node_pos = np.insert(node_pos,np.array(rem_elec_idx)-1,[rem_node[0][0]]*len(rem_elec_idx),axis=0)
        
    #add nodes to mesh
    mesh.add_e_nodes(node_pos-1)#in python indexing starts at 0, in gmsh it starts at 1 
    
    return mesh

#%% column mesh 
def prism_mesh(elec_x,elec_y,elec_z, 
               file_path='column_mesh.geo',
               keep_files=True,
               show_output=True, 
               path='exe', dump=print,
               **kwargs):
    """Make a prism mesh 
    Parameters
    ------------
    elec_x: array like
        electrode x coordinates 
    elec_y: array like 
        electrode y coordinates 
    elec_z: array like 
        electrode z coordinates 
    poly: list, tuple, optional 
        Describes polygon where the argument is 2 by 1 tuple/list. Each entry is the polygon 
        x and y coordinates, ie (poly_x, poly_y)
    z_lim: list, tuple, optional 
        top and bottom z coordinate of column, in the form (min(z),max(z))
    radius: float, optional 
        radius of column
    file_path: string, optional 
        name of the generated gmsh file (can include file path also) (optional)
    cl: float, optional
        characteristic length (optional), essentially describes how big the nodes 
        assocaited elements will be. Usually no bigger than 5. If set as -1 (default)
        a characteristic length 1/4 the minimum electrode spacing is computed.
    elemz: int, optional
        Number of layers in between each electrode inside the column mesh. 
    """
    #error checks 
    if len(elec_x) != len(elec_y):
        raise ValueError("mismatch in electrode x and y vector length, they should be equal")
        
    #make prism mesh command script
    file_name="prism_mesh"
    gw.column_mesh([elec_x,elec_y,elec_z], file_path=file_name, **kwargs)
    #note that this column mesh function doesnt return the electrode nodes 
    
    
    #check directories 
    if path == "exe":
        ewd = os.path.join(
                os.path.dirname(os.path.realpath(__file__)),
                path)
    else:
        ewd = path # points to the location of the .exe 
        # else its assumed a custom directory has been given to the gmsh.exe 
    # handling gmsh
    if platform.system() == "Windows":#command line input will vary slighty by system 
#        cmd_line = ewd+'\gmsh.exe '+file_name+'.geo -3'
        cmd_line = os.path.join(ewd,'gmsh.exe')+' '+file_name+'.geo -3'
    elif platform.system() == 'Darwin':
            winePath = []
            wine_path = Popen(['which', 'wine'], stdout=PIPE, shell=False, universal_newlines=True)#.communicate()[0]
            for stdout_line in iter(wine_path.stdout.readline, ''):
                winePath.append(stdout_line)
            if winePath != []:
                cmd_line = ['%s' % (winePath[0].strip('\n')), ewd+'/gmsh.exe', file_name+'.geo', '-3']
            else:
                cmd_line = ['/usr/local/bin/wine', ewd+'/gmsh.exe', file_name+'.geo', '-3']
    else:
        if os.path.isfile(os.path.join(ewd,'gmsh_linux')): # if linux gmsh is present
            cmd_line = [ewd+'/gmsh_linux', file_name+'.geo', '-3']
        else: # fallback on wine
            cmd_line = ['wine',ewd+'/gmsh.exe', file_name+'.geo', '-3']
        
    if show_output: 
        try:
            p = Popen(cmd_line, stdout=PIPE, shell=False)#run gmsh with ouput displayed in console
        except: # hotfix to deal with failing commits on gitlab's server. 
            cmd_line = ['wine',ewd+'/gmsh.exe', file_name+'.geo', '-3'] # use .exe through wine instead
            p = Popen(cmd_line, stdout=PIPE, shell=False)
        while p.poll() is None:
            line = p.stdout.readline().rstrip()
            dump(line.decode('utf-8'))
    else:
        p = Popen(cmd_line, stdout=PIPE, stderr=PIPE, shell=False)
        p.communicate() # wait to finish

        
    #convert into mesh.dat
    mesh_info = gw.msh_parse(file_path = file_name+'.msh') # read in 3D mesh file
   
    # merge fine with coarse regions
    regions = np.array(mesh_info['parameters'])
    for reg in np.unique(regions)[1:]:
        ie = regions == reg
        regions[ie] = reg - 1
        
    mesh = Mesh(mesh_info['node_x'], # convert output of parser into an object
                mesh_info['node_y'],
                mesh_info['node_z'],
                mesh_info['node_id'],
                mesh_info['elm_id'],
                mesh_info['node_data'],
                mesh_info['cell_type'],
                regions, # this is used as cell_attributes
                mesh_info['parameter_title'],
                mesh_info['original_file_path'])
    
    mesh.moveElecNodes(elec_x,elec_y,elec_z)
    
    if keep_files is False: 
        os.remove(file_name+".geo");os.remove(file_name+".msh")
        
    return mesh 
    
#%% import a custom mesh, you must know the node positions 
def custom_mesh_import(file_path, node_pos=None, order_nodes=True):
    """ 
    Import user defined mesh, currently supports .msh, .vtk and .dat (native to R2/3t)
    format for quad, triangular and tetrahedral meshes. The type of file is guessed from the 
    extension given to the code. 
    
    Parameters
    ---------- 
    file_path : string
        Path to file.
    node_pos : array like, optional
        Array of ints referencing the electrode nodes. If left as none no electrodes 
        will be added to the mesh class. Consider using mesh.move_elec_nodes()
        to add nodes to mesh using their xyz coordinates.
    order_nodes : bool, optional
        Order nodes when importing a mesh  
        
    Returns
    -------
    mesh : class
        mesh class used in ResIPy
        
    """
    if not isinstance(file_path,str):
        raise TypeError("Expected string type argument for 'file_path'")
    
    path,ext = os.path.splitext(file_path)
    if ext == '.vtk':
        mesh = vtk_import(file_path, order_nodes=order_nodes)
    elif ext == '.msh':
        mesh_dict = gw.msh_parse(file_path)

        mesh = Mesh(node_x = mesh_dict['node_x'],
                    node_y = mesh_dict['node_y'],
                    node_z = mesh_dict['node_z'],
                    node_id = mesh_dict['node_id'],
                    elm_id = mesh_dict['elm_id'],
                    node_data = mesh_dict['node_data'],
                    cell_type = mesh_dict['cell_type'],
                    cell_attributes = mesh_dict['parameters'],
                    atribute_title = mesh_dict['parameter_title'],
                    order_nodes = order_nodes)
        
    elif ext == '.dat':
        mesh = dat_import(file_path, order_nodes=order_nodes)   
    elif ext == '.node':
        mesh = tetgen_import(file_path, order_nodes=order_nodes)
    else:
        avail_ext = ['.vtk','.msh','.dat','.node']
        raise ImportError("Unrecognised file extension, available extensions are "+str(avail_ext))
    
    if node_pos is not None:
        mesh.add_e_nodes(np.array(node_pos, dtype=int)) # add electrode nodes to mesh provided by the user
    
    return mesh

#%% ram amount check and is wine installed?. 
#Now for complicated meshes we need alot more RAM. the below function is a os agnostic
#function for returning the amount of total ram. 
# we also need to check wine is installed if running on macOs or linux. 
def systemCheck(show=False):
    """Performs a simple diagnostic of the system, no input commands needed. System
    info is printed to screen, number of CPUs, memory and OS. This check is 
    useful for parallel processing. 
    
    Parameters
    ----------
    show : bool, optional
        If `True`, system specs will be printed.
    
    Returns
    -------
    system_info: dict
        Dictionary keys refer information about the system 
    """
    if show:
        def dump(x):
            print(x)
    else:
        def dump(x):
            pass
    dump("________________System-Check__________________")
    
    totalMemory = '' # incase system can't figure it out!
    num_threads = ''
    OpSys = ''
    #display processor info
    dump("Processor info: %s"%platform.processor())
    num_threads = multiprocessing.cpu_count()
    dump("Number of logical CPUs: %i"%num_threads)
    #this message will display if wine is not installed / detected
    helpful_msg ="""   
This version of ResIPy requires wine to run R2.exe, please consider installing
'wine is not an emulator' package @ https://www.winehq.org/. On linux wine can be found on
most reprositories (ubuntu/debian users can use "sudo apt install wine"). Wine acts as
a compatiblity layer between unix like OS systems (ie macOS and linux) and windows programs. 
    """
    msg_flag = False
    #check operating system 
    OpSys=platform.system()    
    if OpSys=='Darwin':
        dump("Kernel type: macOS")
    else:
        dump("Kernel type: %s"%OpSys)
    #check the amount of ram 
    if OpSys=="Linux":
        p = Popen('free -m', stdout=PIPE, shell=True)
        totalMemory = p.stdout.readlines()[1].split()[1]
        #detect wine 
        p = Popen("wine --version", stdout=PIPE, shell=True)
        is_wine = str(p.stdout.readline())#[0].split()[0]
        if is_wine.find("wine") == -1:
            warnings.warn("Wine is not installed!", Warning)
            msg_flag = True
        else:
            wine_version = is_wine.split()[0].split('-')[1]
            dump("Wine version = "+wine_version)
                          
    elif OpSys=="Windows":
        p = Popen('wmic MEMORYCHIP get Capacity', stdout=PIPE)
        info = p.stdout.readlines()#first entry is the header, subsiquent entries 
        #correspond to dimm slot capacity in bytes 
        totalMemory = 0 # memory returned in binary bytes 
        for i in range(1,len(info)):
            try:
                mem=int(info[i].strip())
                totalMemory += mem
            except ValueError:
                break
        totalMemory = totalMemory/1048576
                
    elif OpSys=='Darwin':
        sysinfo = []
        info = Popen(['system_profiler','SPHardwareDataType'], shell = False, stdout=PIPE, universal_newlines=True)
        for stdout_line in iter(info.stdout.readline, ''):
            sysinfo.append(stdout_line)
        memoryLine = [s for s in sysinfo if any(xs in s for xs in ['Memory'])] 
        totalMemory = re.findall('\\d+', memoryLine[0]) 
        totalMemory = int(totalMemory[0])*1000
        #detect wine
        try: 
            winePath = []
            wine_path = Popen(['which', 'wine'], stdout=PIPE, shell=False, universal_newlines=True)#.communicate()[0]
            for stdout_line in iter(wine_path.stdout.readline, ''):
                winePath.append(stdout_line)
            if winePath != []:
                is_wine = Popen(['%s' % (winePath[0].strip('\n')), '--version'], stdout=PIPE, shell = False, universal_newlines=True)
            else:
                is_wine = Popen(['/usr/local/bin/wine','--version'], stdout=PIPE, shell = False, universal_newlines=True)
            wineVersion = []
            for stdout_line in iter(is_wine.stdout.readline, ""):
                wineVersion.append(stdout_line)
            wine_version = stdout_line.split()[0].split('-')[1]
            dump("Wine version = "+wine_version)
        except:
            warnings.warn("Wine is not installed!", Warning)
            msg_flag = True
        
    else:
        raise OSError("unrecognised/unsupported operating system")
     
    if totalMemory != '':
        totalMemory = int(totalMemory)
        dump("Total RAM available: %i Mb"%totalMemory)
        
        #print some warnings incase the user has a low end PC
        if totalMemory <= 4000:
            warnings.warn("The amount of RAM currently installed is low (<4Gb), complicated ERT problems may incur memory access voilations", Warning)
    
    if num_threads!= '':
        if num_threads <=2:
            warnings.warn("Only one or two CPUs detected, multithreaded workflows will not perform well.", Warning)
            
    if msg_flag:
        dump(helpful_msg)
    
    return {'memory':totalMemory,'core_count':num_threads,'OS':OpSys}