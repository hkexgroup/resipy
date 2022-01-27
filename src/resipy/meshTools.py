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
import os, platform, warnings, psutil
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
from copy import deepcopy

#import R2gui API packages 
import resipy.gmshWrap as gw
from resipy.sliceMesh import sliceMesh # mesh slicing function
import resipy.interpolation as interp

try: 
    from resipy.cext import meshCalc as mc 
except Exception as e:
    print('Could not import meshCalc extension see the following error:')
    print(e)# meshCalc needs to be compiled 
    raise Exception('Could not import meshCalc extension to fix the problem try, '\
                    'updating ResIPy, updating Numpy or recompiling the extension.')

# import pyvista if available
try:
    import pyvista as pv
    try:
        from pyvistaqt import BackgroundPlotter # newer version
    except:
        from pyvista import BackgroundPlotter # older version
    pyvista_installed = True
except:
    pyvista_installed = False
    warnings.warn('pyvista not installed, 3D meshing viewing options will be limited')
    
#%% system status 
# number of cpu cores to use (in functions that can leverage more than 1 core)
global ncores
ncores = 2 
#get gpu info from pyvista 
gpuinfo = None 
# if pyvista_installed:
#     gpuinfo = pv.GPUInfo().renderer


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
    """Determine if a point lies inside a bounding volume.
    
    Parameters
    ----------
    x : array like, float
        x coordinate of query point
    y : array like, float
        y coordinate of query point 
    z : array like, float
        z coordinate of query point
    volume_data : list 
        contains column of poly_data, in the form (polyx, polyy, polyz)
    ray_cast : float, optional
        determines how the far in the x axis a point is ray casted 
    
    Returns
    -------
    inside : boolian, numpy array 
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

#%% check mac version for wine
def getMacOSVersion():
    OpSys=platform.system()    
    if OpSys=='Darwin':
        versionList = platform.mac_ver()[0].split('.')
        macVersion = float(versionList[0] + '.' + versionList[1]) # not getting patch version so xx.xx only
        if macVersion >= 10.15:
            return True
        else:
            return False
        
#%% create mesh object
class Mesh:
    """Mesh class.
    
    Parameters
    ----------
    node_x : list, 1d numpy array
        x coordinates of nodes 
    node_y : list, 1d numpy array
        coordinates of nodes
    node_z : list, 1d numpy array
        z coordinates of nodes 
    connection_matrix: M by N numpy array 
        nodes of element vertices in the form [[node1],[node2],[node3],...], each
        node id should be an integer type. 
    cell_type : list of ints
        code referencing cell geometry (e.g. triangle) according to vtk format
    original_file_path : string, optional
        file path to where the mesh file was originally imported
    regions : optional
        element indexes for a material in the mesh (needs further explanation)
        
    Returns
    ----------
    Mesh : class
    
    Attritbutes 
    ----------
    numnp: int
        Number of nodes 
    numel: int 
        Number of elements 
    node: N by 3 np array 
        Node coordinates (x y z columns)
    connection: N by M np array 
        Connection matrix mapping elements on to thier respective vertices or 
        nodes. Where N is the number of elements and M is the number of element
        vertices (i.e. M = 3 for a triangle mesh)
    elmCentre: N by 3 np array
        Coordinates for the geometric center for each element (x y z columns)
    df: pandas dataframe 
        Contains attributes assocaited with each element inside the mesh, for 
        example resistivity or parameter number 
    Note: 
    ----------
    Something about the attribute_cache 
    """

    #%% mesh creation
    ptdf = None 
    def __init__(self,#function constructs our mesh object. 
                 node_x,#x coordinates of nodes 
                 node_y,#y coordinates of nodes
                 node_z,#z coordinates of nodes 
                 node_data,#nodes of element vertices
                 cell_type,#according to vtk format
                 original_file_path='N/A',
                 order_nodes=True,# order nodes if True, can be computationally expensive 
                 compute_centre=True, # compute cell centres if true, also expensive for big meshes  
                 check2D=True): #check if y and z columns need swapping if mesh is 2D 
        

        #assign variables to the mesh object 
        self.numnp = len(node_x) # ninum >>> use Andy's naming scheme 
        self.numel = node_data.shape[0] # numel 
        self.node = np.array([node_x,node_y,node_z]).T 
        
        dint = 'int64' # connection matrix needs to be in long format 
        if platform.system() == 'Windows':
            dint = np.int32 # avoid windows quirk where type long is actually a 32 bit integer 
        self.connection = np.asarray(node_data,dtype=dint) #connection matrix
        
        self.cell_type = cell_type # cellType
        self.originalFilePath = original_file_path # originalFilePath 
        self.eNodes  = None # node number that electrodes occupy 
        self.elec = None
        self.iremote = None # specify which electrode is remote
        self.cax = None # store mesh.show() output for faster mesh.draw()
        self.zone = np.ones(self.numel) # by default all in the same zone # remove? 
        self.elmCentre = None  # remove? 
        
        #mesh element attributes 
        df={'param': np.arange(self.numel)+1,
            'elm_id': np.arange(self.numel)+1,
            #'zones':np.ones(self.numel) + 1,
            'region': np.ones(self.numel),
            'cellType':np.ones(self.numel)*cell_type[0]}# store attributes values per cell
        self.df= pd.DataFrame(df) # rename to df 
        
        self.mesh_title = "2D_R2_mesh"
        self.no_attributes = 0
        
        #mesh calculations 
        self.neigh_matrix = None # neighbour matrix, not usually needed unless for 3d tetrahedra problems 
        self.tri_combo = None
        self.NsizeA = None # Nconnec + numnp  
        self.fconm = None # finite element conductance matrix
        
        #decide on number of dimensions
        if max(node_y) - min(node_y) < 1e-16 or max(node_z) - min(node_z) < 1e-16: # mesh is probably 2D 
            self.ndims=2
        else:
            self.ndims=3
            self.mesh_title = '3D_R3t_mesh'
            
        if self.ndims == 2 and max(node_z) - min(node_z) == 0 and check2D: # then the mesh is lightly to be a 2D mesh but with elevation in Y axis
            warnings.warn('Y and Z columns in mesh are flipped so that elevation is in the Z axis, to reverse this use mesh.flipYZ()')
            self.flipYZ()
        
        if compute_centre: #compute the centre of mesh cells
            self.cellCentres()
        if order_nodes: # order nodes 
            self.orderNodes()
            
        self.surfaceMesh = None # surface of mesh 
            
        #point data attributes 
        ptdf = {'x':node_x,'y':node_y,'z':node_z}
        self.ptdf = pd.DataFrame(ptdf)
    
    def copy(self):
        """Return a copy of mesh object. 
        """
        # mesh = deepcopy(self) # not always working! below too not always working!!!
        
        mesh = Mesh(self.node[:,0],
                    self.node[:,1],
                    self.node[:,2],
                    self.connection,
                    self.cell_type,
                    self.originalFilePath,
                    order_nodes=False,# should be already ordered
                    compute_centre=False,# should be already computed 
                    check2D=False)# should be already computed 
        mesh.elmCentre = self.elmCentre
        mesh.df = self.df.copy()
        
        #electrode handling
        try:
            mesh.eNodes = self.eNodes
            if self.eNodes is not None:
                mesh.setElecNode(self.eNodes)
        except:
            pass
        mesh.iremote = self.iremote
        mesh.zone = self.zone
        
        return mesh 
    
    
    def flipYZ(self):
        """ Make Y column Z column and vice versa, this is useful for when 
        dealing with 2D meshes. Dimensions are modified in place. 
        """
        node_y_cache = self.node[:,1].copy()
        node_z_cache = self.node[:,2].copy()
        
        self.node[:,1] = node_z_cache
        self.node[:,2] = node_y_cache
        
        self.cellCentres() # recompute cell centres 
                
    
    #%% mesh attribute handling 
    def setElecNode(self, e_nodes, iremote =None):
        """Assign node numbers to electrodes. 
        
        Parameters
        ------------
        e_nodes: array like
            array of ints which index the electrode nodes in a mesh
        iremote: array like, optional
            Array of bool, if True then indexed electrode is remote. 
        """
        self.eNodes  = e_nodes
        self.elec = self.node[e_nodes,:]
        if iremote is None:
            self.iremote = np.array([False]*self.elec.shape[0],dtype=bool)
        else:
            self.iremote = iremote 
            
    def setElec(self,elec_x,elec_y,elec_z):
        if len(elec_x) != len(elec_y) or len(elec_x) != len(elec_z):
            raise ValueError('Mismatch in the electrode array lengths when setting the electrodes in mesh class')
        self.elec = np.array([elec_x,elec_y,elec_z]).T
        if self.eNodes is not None:
            self.moveElecNodes(elec_x,elec_y,elec_z,debug=False)
        self.iremote = np.array([False]*self.elec.shape[0],dtype=bool)
    
    def add_e_nodes(self, e_nodes):
        warnings.warn('add_e_nodes is deprecaited: use setElecNode instead', DeprecationWarning)
        self.setElecNode(e_nodes)
    
    #add some functions to allow adding some extra attributes to mesh 
    def addSensitivity(self,values):#sensitivity of the mesh
        if len(values)!=self.numel:
            raise ValueError("The length of the new attributes array does not match the number of elements in the mesh")
        self.sensitivities = values
        
    def file_path(self):
        """Returns the file path from where the mesh was imported
        """
        return(format(self.originalFilePath))
       
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
        
    def type2FaceNo(self):
        """Converts vtk cell types into number of faces each element has
        """
        if int(self.cell_type[0])==5:#then elements are triangles
            return 1
        elif int(self.cell_type[0])==8 or int(self.cell_type[0])==9:#elements are quads
            return 1
        elif int(self.cell_type[0]) == 11: # elements are voxels
            return 8
        elif int(self.cell_type[0]) == 10:# elements are tetrahedra 
            return 4
        elif int(self.cell_type[0]) == 13: # elements are 3d wedges 
            return 5
        #add element types as neccessary 
        else:
            print("WARNING: unrecognised cell type")
            return 0
        
    def summary(self,flag=True):
        """Prints summary information about the mesh
        """
        self.no_attributes = len(self.df.keys())
        #returns summary information about the mesh, flagto print info, change to return string
        out = "\n_______mesh summary_______\n"
        out += "Number of elements: %i\n"%int(self.numel)
        out += "Number of nodes: %i\n"%int(self.numnp)
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
        return self.summary(flag=False) + self.showAvailAttr(flag=False)
            
    def add_attribute(self,values,key):
        warnings.warn('add_attribute is depreciated, use addAttribute instead', DeprecationWarning)
        self.addAttribute(values, key)
        
    def addAttribute(self,values,key):
        """Add a new/update attribute to mesh. 
        
        Parameters
        ------------
        values: array like
            must have a length which matches the number of elements. Discrete 
            values which map to the elements. 
        key: str
            Name of the attribute, this will be used to reference to the values
            in other mesh functions. 
        """
        if len(values)!=self.numel:
            raise ValueError("The length of the new attributes array (%i) does not match the number of elements in the mesh (%i)"%(len(values),self.numel))
        self.no_attributes += 1
        self.df[key]=values #allows us to add an attributes to each element.
        
    def addPtAttribute(self,values,key):
        """Associate attributes with the mesh nodes  
        
        Parameters
        ------------
        values: array like
            Must have a length which matches the number of nodes. Discrete 
            values which map to the nodes. 
        key: str
            Name of the attribute, this will be used to reference to the values
            in other mesh functions. 
        """
        if len(values)!=self.numnp:
            raise ValueError("The length of the new attributes array (%i) does not match the number of nodes in the mesh (%i)"%(len(values),self.numnp))
        self.ptdf[key]=values #allows us to add an attributes to each element.
        
    def show_avail_attr(self,flag=True):
        warnings.warn("show_avail_attr is depreciated, use showAvailAttr instead", DeprecationWarning)
        self.showAvailAttr(flag)
        
    def showAvailAttr(self,flag=True):
        """Show available attributes in mesh.df. 
        """
        out = '\n______cell attributes_____\n'
        try: 
            for i,key in enumerate(self.df.keys()):
                out += key + '\n'
        except:
            out += "None\n"
        if flag:
            print(out)
        else:
            return out
    
    #%% mesh calculations 
    def orderNodes(self, return_count=False):
        """Order mesh nodes in clockwise fashion 
        
        Parameters
        -----------
        return_count:bool, optional
            Function returns the number of reordered elements, default is False. 
        """
        con_mat = self.connection
        dint = self.connection.dtype
        con_mata = self.connection.copy() # new object of the connection matrix to reference
        node_x = self.node[:,0]
        node_y = self.node[:,1]
        node_z = self.node[:,2]
        count = 0
        
        if int(self.cell_type[0])==5:#then elements are triangles, now vectorised
            p = np.array((node_x[con_mata[:,0]], node_z[con_mata[:,0]])).T
            q = np.array((node_x[con_mata[:,1]], node_z[con_mata[:,1]])).T
            r = np.array((node_x[con_mata[:,2]], node_z[con_mata[:,2]])).T
            ccw = interp.ccwv(p, q, r)
            ds = ccw==1 # do switch 
            ns = ccw!=1 # no switch
            con_mat = np.zeros_like(con_mata)
            con_mat[:,0] = con_mata[:,0]
            con_mat[ds,1] = con_mata[ds,2]
            con_mat[ns,1] = con_mata[ns,1]
            con_mat[ds,2] = con_mata[ds,1]
            con_mat[ns,2] = con_mata[ns,2]
            count = np.sum(ds)
                    
        elif int(self.cell_type[0])==8 or int(self.cell_type[0])==9:#elements are quads
            con_mat, count = mc.orderQuad(self.connection,
                                          self.node) # call cython 
                    
        elif int(self.cell_type[0]) == 11: # elements are voxels
            #Node ordering scheme not avialable with this mesh type
            return
        
        elif int(self.cell_type[0]) == 10:# elements are tetrahedra 
            con_mat, count, _ = mc.orderTetra(self.connection,
                                              self.node,
                                              ncores) # call cython 
                    
        elif int(self.cell_type[0]) == 13: # elements are 3d wedges 
            for i in range(self.numel):
                n1=(node_x[con_mat[i][0]],node_y[con_mat[i][0]],node_z[con_mat[i][0]])#define node coordinates
                n2=(node_x[con_mat[i][1]],node_y[con_mat[i][1]],node_z[con_mat[i][1]])
                n3=(node_x[con_mat[i][2]],node_y[con_mat[i][2]],node_z[con_mat[i][2]])

                #see if top of triangle is counter-clockwise
                if interp.ccw(n1,n2,n3) == 1: #points are clockwise and therefore need swapping round
                    count += 1
                    con_mat[i][1] = con_mata[i][0] # correct the top of the triangle 
                    con_mat[i][0] = con_mata[i][1]
                    con_mat[i][3] = con_mata[i][4] # correct the bottom of the triangle 
                    con_mat[i][4] = con_mata[i][3]

        
        self.connection = np.asarray(con_mat,dtype=dint)
        
        if return_count:
            return count
        
    def orderElem(self, param=None):
        """Order elements based on the parameter number. Ideally parameters 
        should be concurrent, and fixed elements should go at the base of the 
        connection matrix
        """
        if param is None:
            param = self.df['param'].values
        else:
            if len(param)!= self.numel:
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
           
        con_mata = self.connection.copy()
        con_mat_fix = con_mata[order,:]
        paramFixed = np.array(param)[order]
            
        for key in self.df.keys():
            self.df[key] = self.df[key].values[order]
        
        self.connection = con_mat_fix
        self.cellCentres() # recompute them too
        self.addAttribute(paramFixed,'param')
        
        
    def resetParam(self):
        """Reorder parameters into consective ascending order 
        """
        param = np.arange(self.numel) + 1
        self.addAttribute(param,'param')
        
    
    def cellCentres(self):
        """A numpy-based approximation of cell centres for 2D and 3D elements. 
        It's calculated from the mean of cell x y z node coordinates 
        """
        #print('Calculating cell centres')
        con_mat = self.connection
        vertx = self.node[:,0][con_mat]
        verty = self.node[:,1][con_mat]
        vertz = self.node[:,2][con_mat]

        elm_centre = np.zeros((self.numel,3))
        
        elm_centre[:,0] = np.mean(vertx,axis=1)
        elm_centre[:,1] = np.mean(verty,axis=1)
        elm_centre[:,2] = np.mean(vertz,axis=1)
        self.elmCentre = elm_centre
        self.df['X'] = elm_centre[:,0]
        self.df['Y'] = elm_centre[:,1]
        self.df['Z'] = elm_centre[:,2]
    
    
    def cellArea(self):
        """Compute the element areas, or in the case of 3D meshes compute the 
        cell volumes. 
        """
        con_mat = self.connection
        elm_area=[0]*self.numel
        node_x = self.node[:,0]
        node_y = self.node[:,1]
        node_z = self.node[:,2]

        if int(self.cell_type[0])==5:#then elements are triangles
            for i in range(self.numel):
                n1=(node_x[con_mat[i][0]], node_z[con_mat[i][0]])#define node coordinates
                n2=(node_x[con_mat[i][1]], node_z[con_mat[i][1]])
                n3=(node_x[con_mat[i][2]], node_z[con_mat[i][2]])
                #compute area (for a triangle this is 0.5*base*height)
                base=(((n1[0]-n2[0])**2) + ((n1[1]-n2[1])**2))**0.5
                mid_pt=((n1[0]+n2[0])/2,(n1[1]+n2[1])/2)
                height=(((mid_pt[0]-n3[0])**2) + ((mid_pt[1]-n3[1])**2))**0.5
                elm_area[i] = 0.5*base*height
                
        elif int(self.cell_type[0])==8 or int(self.cell_type[0])==9:#elements are quads
            for i in range(self.numel):
                p = self.node[con_mat[i]].T
                dx = abs(max(p[0]) - min(p[0]))
                dz = abs(max(p[1]) - min(p[1]))
                elm_area[i] = dx*dz
                
        elif int(self.cell_type[0]) == 11: # elements are voxels
            for i in range(self.numel):
                p = self.node[con_mat[i]].T
                #compute volume (which is a bit of an approximation)
                dx = abs(max(p[0]) - min(p[0]))
                dy = abs(max(p[1]) - min(p[1]))
                dz = abs(max(p[2]) - min(p[2]))
                elm_area[i] = dx*dy*dz

        elif int(self.cell_type[0]) == 10:# elements are tetrahedra 
            warnings.warn('Area calculation for tetrahedral meshes is not yet validated')
            for i in range(self.numel):
                p = self.node[con_mat[i]].T
                #find apex of tetra (solve via dot product)
                P = p[:,0] # point 1 
                Q = p[:,1] # point 2 
                R = p[:,2] # point 3
                pq = Q - P # p to q vector 
                pr = R - P # p to r vector 
                v = np.cross(pq, pr) # cross product 
                A = v[0] # equation of plane parameters 
                B = v[1]
                C = v[2]
                D = -(A*P[0] + B*P[1] + C*P[2])
                S = p[:,3] # point 4
                #height of tetra 
                height = np.abs(A*S[0] + B*S[1] + C*S[2] + D) / np.sqrt(A**2 + B**2 + C**2)
                #area of base triangle 
                PQ = np.linalg.norm(pq) # magnitude / norm of vectors 
                PR = np.linalg.norm(pr)
                theta = np.arccos(np.dot(pq,pr) / (PQ*PR)) # find angle between vectors 
                area = 0.5 * PQ * PR * np.sin(theta) #area of baseline triangle 
                
                #according to http://mathcentral.uregina.ca/QQ/database/QQ.09.03/peter2.html
                #the volume is calculated as 1/3 * height * area of triangle 
                elm_area[i] = (1/3)*area*height
                
        elif int(self.cell_type[0]) == 13: # elements are 3d wedges 
            for i in range(self.numel):
                n1=(node_x[con_mat[i][0]], node_y[con_mat[i][0]], node_z[con_mat[i][0]])#define node coordinates
                n2=(node_x[con_mat[i][1]], node_y[con_mat[i][1]], node_z[con_mat[i][1]])
                n3=(node_x[con_mat[i][2]], node_y[con_mat[i][2]], node_z[con_mat[i][2]])
                n4=(node_x[con_mat[i][3]], node_y[con_mat[i][3]], node_z[con_mat[i][3]])
                n5=(node_x[con_mat[i][4]], node_y[con_mat[i][4]], node_z[con_mat[i][4]])
                n6=(node_x[con_mat[i][5]], node_y[con_mat[i][5]], node_z[con_mat[i][5]])
                #compute wedge volume by computing face area first
                base=(((n1[0]-n2[0])**2) + ((n1[1]-n2[1])**2))**0.5
                mid_pt=((n1[0]+n2[0])/2,(n1[1]+n2[1])/2)
                height=(((mid_pt[0]-n3[0])**2) + ((mid_pt[1]-n3[1])**2))**0.5
                area = 0.5*base*height
                p = np.array((n1,n2,n3,n4,n5,n6)).T
                dz = abs(max(p[2]) - min(p[2]))
                elm_area[i] = area * dz
                
        if self.ndims == 2:
            self.df['Area'] = elm_area
        else:                
            self.df['Volume'] = elm_area

            
    def computeNeigh(self): # fix me 
        """Compute element neighbour matrix
        """            
        if self.ndims == 2: #2d mesh 
            self.neigh_matrix, self.tri_combo = mc.neigh2d(self.connection,1)
        elif self.ndims == 3: #3d mesh 
            if self.type2VertsNo() == 4:# tetra mesh 
                self.neigh_matrix, self.tri_combo = mc.neigh3d(self.connection,1,ncores)
            elif self.type2VertsNo() == 6:  # prism mesh 
                self.neigh_matrix, self.tri_combo = mc.neighPrism(self.connection,1,ncores)


    def computeNconnec(self):
        cell_type = self.cell_type[0]
        self.NsizeA, self.fconm = mc.conductanceCall(self.connection, self.numnp,
                                                     cell_type,ncores)
        
    def refine(self):
        """Refine the mesh into smaller elements 
        """
        if self.ndims==2 and self.type2VertsNo()==3:
            #then its a triangle mesh
            return self.splitTri()
        elif self.ndims==2 and self.type2VertsNo()==4:
            # then its a quad mesh (returns a new mesh)
            return self.quad2tri()
        elif self.ndims==3 and self.type2VertsNo()==4:
            #then its a tetra mesh 
            return self.splitTetra()
        else:
            print('Sorry not implimented for this mesh type yet')
            return 
        
        
    def splitTri(self, param=None): # fix me 
        """Refine triangles by splitting them into 4 smaller triangles 
        """
        #error checks 
        if self.ndims==3:
            raise ValueError("This kind of mesh splitting isn't available for 3D meshes")
        elif self.type2VertsNo() == 4: #not a triangle mesh 
            raise Exception("Sorry mesh splitting not avialable for this mesh type")
        
        #see if parameter already assigned 
        if param is None:
            if 'param' not in self.df.keys():
                self.df['param'] = 1 + np.arange(self.numel)
            param = self.df['param'].copy()
        else:
            if len(param)!= self.numel:
                raise ValueError('The parameter array does not match the number of elements')
                return
                    
        (new_con_mat, new_node,
         nnum_elms, nnum_nodes) = mc.splitTri(self.connection, self.node)
             
        
        nmesh = self.copy() 
        nmesh.connection = new_con_mat
        nmesh.node = new_node
        nmesh.numel = nnum_elms
        nmesh.numnp = nnum_nodes
        
        if self.df is not None:
            df = self.df.copy()
            new_df = {}
            for key in self.df.keys():
                a = np.repeat(df[key].values,4)
                new_df[key] = a
                
        if self.ptdf is not None:
            new_ptdf = {}
            new_ptdf['x'] = new_node[:,0]
            new_ptdf['y'] = new_node[:,1]
            new_ptdf['z'] = new_node[:,2]
            nmesh.ptdf = pd.DataFrame(new_ptdf)
            
        nmesh.df = pd.DataFrame(new_df)
        nmesh.cellCentres()
        nmesh.orderNodes() 
        return nmesh 
             
        
    def splitTetra(self, param=None): # fix me 
        """Refine tetrahedra by splitting them in six smaller tetrahedra.
        """
        #error checks 
        if self.ndims == 2:
            raise ValueError("This kind of mesh splitting isn't available for 2D meshes")
        if self.type2VertsNo() == 6:#not a tetrahedra 3d mesh 
            raise Exception("This kind of mesh splitting isn't available for meshes of type 'prism'")
            
        #see if parameter already assigned 
        if param is None:
            if 'param' not in self.df.keys():
                self.df['param'] = 1 + np.arange(self.numel)
            param = self.df['param'].copy()
        else:
            if len(param)!= self.numel:
                raise ValueError('The parameter array does not match the number of elements')
                return
        
        (new_con_mat, new_node,
         nnum_elms, nnum_nodes) = mc.splitTetra(self.connection, self.node)
         
        nmesh = self.copy() 
        nmesh.connection = new_con_mat
        nmesh.node = new_node
        nmesh.numel = nnum_elms
        nmesh.numnp = nnum_nodes
        #nmesh.cell_attributes = np.repeat(self.cell_attributes,8)
        
        if self.df is not None:
            df = self.df.copy()
            new_df = {}
            for key in self.df.keys():
                a = np.repeat(df[key].values,8)
                new_df[key] = a
                
        nmesh.df = pd.DataFrame(new_df)
                
        if self.ptdf is not None:
            new_ptdf = {}
            new_ptdf['x'] = new_node[:,0]
            new_ptdf['y'] = new_node[:,1]
            new_ptdf['z'] = new_node[:,2]
            nmesh.ptdf = pd.DataFrame(new_ptdf)
        
        nmesh.cellCentres()
        nmesh.orderNodes()
        
        #mesh calculations - ensure all these are set to none so they are recalculated 
        nmesh.neigh_matrix = None # neighbour matrix, not usually needed unless for 3d tetrahedra problems 
        nmesh.tri_combo = None
        nmesh.NsizeA = None # Nconnec + numnp  
        nmesh.fconm = None # finite element conductance matrix
        return nmesh
    
    def quad2tri(self):
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
        
        con_mata = self.connection.copy() # get connection matrix 
        con1 = con_mata[:,[0,1,2]] # retrieve vertices of columns 
        con2 = con_mata[:,[0,3,2]] # and split to make right angle triangles 
        out_mat = np.zeros((self.numel*2,3),dtype=int)       
        out_mat[0::2,:] = con1
        out_mat[1::2,:] = con2 
        
        tri_mesh = self.copy()
        tri_mesh.connection = out_mat
        tri_mesh.numel = self.numel*2
        tri_mesh.cell_type=[5]
        
        tri_df = {}
        
        if self.df is not None:
            calculated = ['X','Y','Z']
            for key in self.df.keys():
                if key not in calculated:
                    a = np.repeat(self.df[key],2) # repeat cell attributes 
                    tri_df[key] = a
        
        tri_mesh.df = pd.DataFrame(tri_df) #reassign dataframe 
        tri_mesh.orderNodes() # order the nodes 
        tri_mesh.cellCentres() # calculate cell centres 
        return tri_mesh
        
        
    def elemDist(self):
        """Work out the distance of each cell in the mesh to its closest electrode
        """
        try: 
            elec = self.elec.copy()
        except:
            warnings.warn('No electrodes have been set in mesh class, aborting ... ')
            return
        points = self.elmCentre.copy()
        
        tree = cKDTree(elec) ### setup tree
        dist,idx = tree.query(points) ### >> maps to points to nearest electrode
        self.addAttribute(dist,'cell_distance')
        return dist
    
    def extractSurface(self, return_idx =False, post_neigh_check=True): 
        """ Extract the surface of a triangle or tetrahedral mesh. Ouput of 
        function will depend on mesh type. 
        
        Parameters
        -----------
        return_idx: bool
            Return the indices of elements on the top of the mesh (default is false)
        post_neigh_check: bool 
            Perform a post processing step to check elements have neighbours, 
            this is useful for meshes with uneven sides where false postives 
            can crop up. 
    
        Returns
        -------
        mesh: class
            2d faces of the top of the mesh, if the input mesh is a tetrahedral 
            mesh
        (x,z): tuple
            1D faces of the top of the mesh, if the input is a triangular mesh 
        idx: array like 
            if return_idx == True an array of int
        """
        if self.ndims == 2:
            typ = 2
            if self.type2VertsNo() == 4: 
                qmesh = self.quad2tri() # get the triangle mesh instead 
                return qmesh.extractSurface() # run the same function 

        elif self.ndims == 3: 
            typ = 3
            if self.type2VertsNo() == 6 or self.type2VertsNo() == 8:#not a tetrahedra 3d mesh 
                raise Exception("Sorry surface extraction is not available for this type of mesh")
        
        con_mat = self.connection
        if self.neigh_matrix is None: # compute neighbour matrix 
            self.computeNeigh() # this will find the element neighbours and the elements which lie on the outside of the mesh! 
        neigh = self.neigh_matrix
        
        out_elem = np.min(neigh, axis=1) == -1 # elements which have a face on the outside of the mesh 
        neigh_trunc = neigh[out_elem]
        con_trunc = con_mat[out_elem]           
        
        node_x = self.node[:,0]
        node_y = self.node[:,1]
        node_z = self.node[:,2]
        
        if typ == 2: ### Extract 2D faces of top of the mesh ###
            xm = self.elmCentre[:,0]
            zm = self.elmCentre[:,2]
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
            if return_idx:
                return (uni_x,uni_z), attr_idx
            else:
                return (uni_x,uni_z)
            
        else:  ### Extract 3D faces of top of the mesh ###
            fcon, outelem = mc.faces3d(self.connection,
                                       neigh) # first get external faces 
            # outelem = np.array(outelem)
            points = self.elmCentre[outelem]#get middle points of elements with external face 
                
            # extract surface cells             
            ochecka = mc.surfaceCall(fcon, self.node, points, num_threads=ncores) # send to c extension 
        
            ikeep = ochecka == 1
            nmesh = Mesh(node_x, # make new mesh 
                         node_y, 
                         node_z, 
                         node_data = fcon[ikeep,:], 
                         cell_type = [5], 
                         order_nodes=False,
                         compute_centre=False,
                         check2D=False)
            
            #sort out cell attributes
            df = self.df.copy()
            d = {}
            for key in df.keys():
                a = df[key].values
                d[key] = a[outelem[ikeep]]
                
            nmesh.df = pd.DataFrame(d)
            nmesh.cellCentres()
            
            if post_neigh_check: 
                # remove faces without 3 nieghbours 
                nneigh = mc.neigh2d(nmesh.connection,0)
                ikeep2 = np.min(nneigh,axis=1) > -1
                nmesh = nmesh.filterIdx(ikeep2)
                
                # remove faces with 1 or less neighbours 
                nneigh = mc.neigh2d(nmesh.connection,0)
                ikeep3 = np.count_nonzero(nneigh+1,axis=1) > 1
                nmesh = nmesh.filterIdx(ikeep3) 
              
            nmesh.__rmexcessNodes() # remove excess nodes which are not used 
            
            if return_idx:
                return nmesh, outelem[ikeep] 
            # might cuase a bug where post_neigh_check == True
            else:
                return nmesh
        
    #%% Truncating the mesh 
    def __rmexcessNodes(self):
        """ Remove any nodes are not inside the connection matrix
        """
        con_mat = self.connection.copy()
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
        node = self.node.copy()
        self.node = node[uni_nodes,:]
        new_conf = new_nodes[map_idx]
        new_con = new_conf.reshape(shp)
        self.connection = new_con
            
        self.numnp = len(uni_nodes)
        
        #sort point dataframe 
        if self.ptdf is not None:
            ptdf = self.ptdf.copy()   
            a = np.array([False]*len(ptdf),dtype=bool)
            a[uni_nodes] = True
            new_ptdf = ptdf[a].reset_index()
            self.ptdf = new_ptdf 
        
        
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
        centroids = np.c_[self.elmCentre[:,0], self.elmCentre[:,2]]
        i2keep = path.contains_points(centroids) # this is the fastest way to check if point are in polygon 
        
        # filter element-based attribute
        ec = self.elmCentre.copy()[i2keep,:]
        
        nmesh = deepcopy(self)#.copy()
        nmesh.elmCentre = ec 
        nmesh.connection = self.connection[i2keep,:]
        nmesh.elm_id = np.arange(1,self.numel+1)[i2keep]
        #nmesh.cell_attributes = np.array(self.cell_attributes)[i2keep]
        
        df = self.df.copy()
        nmesh.df = df[i2keep].reset_index()
        del nmesh.df['index']
        nmesh.numel = np.sum(i2keep)
        
        nmesh.__rmexcessNodes()         
        return nmesh
        
    # return a truncated mesh (for 3D)
    def truncateMesh(self, xlim=None, ylim=None, zlim=None):
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
        nmesh: Class
            New instance of mesh class which is truncated 
        """
        if self.elec is not None:
            if xlim is None:
                xlim=[min(self.elec[:,0]), max(self.elec[:,0])]
            if ylim is None:
                ylim=[min(self.elec[:,1]), max(self.elec[:,1])]
        else:
            if xlim is None:
                xlim=[min(self.node[:,0]), max(self.node[:,0])]
            if ylim is None:
                ylim=[min(self.node[:,1]), max(self.node[:,1])]
        if zlim is None:
            zlim=[min(self.node[:,2]), max(self.node[:,2])]
        # why not truncate as well to electrode extent?
        
        # elm_x = self.elmCentre[:,0]
        # elm_y = self.elmCentre[:,1]
        # elm_z = self.elmCentre[:,2]
        # ie = (elm_x > xlim[0]) & (elm_x < xlim[1]) &\
        #      (elm_y > ylim[0]) & (elm_y < ylim[1]) &\
        #      (elm_z > zlim[0]) & (elm_z < zlim[1])
             
        # we keep all the elements which have at least 1 node inside
        # the zone of interest
        ie = (self.node[self.connection,0] > xlim[0]).any(1) &\
             (self.node[self.connection,1] > ylim[0]).any(1) &\
             (self.node[self.connection,2] > zlim[0]).any(1) &\
             (self.node[self.connection,0] < xlim[1]).any(1) &\
             (self.node[self.connection,1] < ylim[1]).any(1) &\
             (self.node[self.connection,2] < zlim[1]).any(1)
                     
        nmesh = self.copy() # make a new mesh object with fewer elements 
        nmesh.df = nmesh.df[ie].reset_index(drop=True)
        nmesh.df['elm_id'] = 1 + np.arange(nmesh.df.shape[0])
        nmesh.numel = nmesh.df.shape[0]
        nmesh.elmCentre = self.elmCentre[ie,:]
        nmesh.connection = nmesh.connection[ie,:]
        
        nmesh.__rmexcessNodes() # remove the excess nodes 
            
        return nmesh # return truncated mesh 
    
    
    def threshold(self, attr=None, vmin=None, vmax=None):
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
        if self.df is None:
            raise EnvironmentError("No 'df' varaible exists for the mesh class!")
        elif attr not in self.df.keys():
            raise ValueError("Specified attribute has not been defined")
            
        X = self.df[attr].values
        
        if vmin is None:
            vmin = np.min(X)
        if vmax is None:
            vmax = np.max(X)
            
        in_elem = (X >= vmin) & (X <= vmax) # index of elements to keep 
        
        temp_con_mat = self.connection.copy() #temporary connection matrix which is just the elements inside the box
        
        #new_attr = np.array(self.cell_attributes)
        
        elm_id = np.arange(self.numel)+1

        #truncate the attribute table down to the inside elements 
        new_df = self.df[in_elem]
        
        nmesh = self.copy() # make a new mesh object with fewer elements 
        
        nmesh.df = new_df
        nmesh.numel = len(elm_id[in_elem])
        # nmesh.cell_attributes = new_attr[in_elem]
        # nmesh.elm_id = elm_id[in_elem]
        nmesh.elmCentre = self.elmCentre[in_elem,:]
        
        new_con_mat = temp_con_mat[in_elem,:]
        nmesh.connection = new_con_mat     
        
        nmesh.__rmexcessNodes() # remove the excess nodes 
            
        return nmesh # return truncated mesh 
    
    def filterIdx(self, in_elem):
        """Filter mesh down on the basis of element number / index
        
        Parameters
        ------------
        in_elem: array like
            array of bool, True means to keep the element 
        attr: string
            Name of attribute to threshold by 
        
        Returns
        ------------
        mesh: Class
            New instance of mesh class which is filtered 
        """
        if self.df is None:
            raise EnvironmentError("No 'df' varaible exists for the mesh class!")
        
        if len(in_elem) != self.numel:
            raise ValueError('Number of filtered indices does not match number of mesh elements')
        
        temp_con_mat = self.connection.copy() #temporary connection matrix which is just the elements inside the box
        
        #new_attr = np.array(self.cell_attributes)
        
        elm_id = np.arange(self.numel)+1

        #truncate the attribute table down to the inside elements 
        new_df = self.df[in_elem].copy() # TODO copy needed? (gb)
        
        nmesh = self.copy() # make a new mesh object with fewer elements 
        
        nmesh.df = new_df
        nmesh.numel = len(elm_id[in_elem])
        # nmesh.cell_attributes = new_attr[in_elem]
        # nmesh.elm_id = elm_id[in_elem]
        nmesh.elmCentre = self.elmCentre[in_elem,:]
        
        new_con_mat = temp_con_mat[in_elem,:]
        nmesh.connection = new_con_mat     
        
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
        
        unix = np.unique(self.node[:,0]) # unique x values in the x node coordinates 
        if topo is None: # find the y column is a little challanging without knowing the original topography 
            uniz = np.unique(self.node[:,2]) # if you dont know any better then just use the unique z values 
        else:
            uniz = topo # ideally use mesh topography 
        e_nodes = self.eNodes 
        colx = [0]*len(e_nodes) # column indexes for x coordinates 
        colz = [0]*len(e_nodes) # column indexes for z coordinates 
        node_x = self.node[:,0]
        node_z = self.node[:,2]
        for i in range(len(e_nodes)):
            x = node_x[e_nodes[i]] # get node x coordinate 
            z = node_z[e_nodes[i]]
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
         
        """
        #assign coordinate arrays 
        x_old = look_up_mesh.elmCentre[:,0]
        x_new = self.elmCentre[:,0]
        y_old = look_up_mesh.elmCentre[:,1]
        y_new = self.elmCentre[:,1]
        z_old = look_up_mesh.elmCentre[:,2]
        z_new = self.elmCentre[:,2]
        i_old = np.ones(look_up_mesh.numel) # dummy parameter 
        #do look up 
        if self.ndims==3:
            i_new, idxes = interp.nearest3d(x_new,y_new,z_new,
                                            x_old,y_old,z_old,i_old,
                                            return_idx=True,
                                            num_threads=ncores)
        elif self.ndims==2:
            i_new, idxes = interp.nearest(x_new,z_new,
                                          x_old,z_old,i_old,
                                          return_idx=True,
                                          num_threads=ncores) 
                    
        #look up values from look up mesh     
        look_up_cache = look_up_mesh.df.copy()
        new_df = {}
        for key in look_up_cache.keys():
            x = look_up_cache[key].values[idxes]
            new_df[key] = x
            
        #keep any attributes already in own dataframe 
        for key in self.df.keys():
            if key not in look_up_cache.keys():
                new_df[key] = self.df[key].values
        
        self.df = pd.DataFrame(new_df) # map indexes using dataframe 
            
    def transMesh(self,x,y,z):
        """Translate mesh by x y z coordinate
        """

        self.node[:,0] += x
        self.node[:,1] += y
        self.node[:,2] += z
        self.cellCentres()
        
    def trans_mesh(self,x,y,z):
        warnings.warn("trans_mesh is depreciated, use transMesh instead", DeprecationWarning)
        self.transMesh(x, y, z)
        
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
        if self.ndims == 3:
            raise ValueError('Cant clip contour on a 3D mesh')
        # mask outer region
        xmin = np.min(self.node[:,0])
        xmax = np.max(self.node[:,0])
        zmin = np.min(self.node[:,2])
        zmax = np.max(self.node[:,2])
        
        (xsurf, zsurf) = self.extractSurface() # extract 2d mesh surface 
        if maxDepth is not None:
            xfmd, zfmd = xsurf[::-1], zsurf[::-1] - maxDepth
            verts = np.c_[np.r_[xmin, xmin, xsurf, xmax, xmax, xfmd, xmin],
                          np.r_[zmin, zmax, zsurf, zmax, zmin, zfmd, zmin]]
        else:
            verts = np.c_[np.r_[xmin, xmin, xsurf, xmax, xmax, xmin],
                          np.r_[zmin, zmax, zsurf, zmax, zmin, zmin]]     
             
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
             darkMode = False,
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
        darkMode : bool, optional
            If True, electrodes will be plotted in white, else black
        
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
        
        if self.iremote is None and self.elec is not None:
            try:
                iremote = np.zeros(self.elec.shape[0], dtype=bool)
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
             darkMode=darkMode,
             **kwargs) # show 3D mesh instead 
            return # exit 2D mesh show function 
            
        
        # decide which attribute to plot, we may decide to have other attritbutes! 
        keys = self.df.keys()
        if attr not in keys and not attr is None:
            attr = None
            warnings.warn('Chosen attribute not found in mesh class')
        if attr is None:
            if 'region' in keys:
                attr = 'region'
            else: 
                attr = keys[0]

        X = np.array(self.df[attr])
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
        
        # if no dimensions are given then set the plot limits to edge of mesh
        try: 
            elec_x = self.elec[:,0][~iremote]
            if xlim==None:
                xlim=[min(elec_x),max(elec_x)]
            if zlim==None:
                doiEstimate = 2/3*np.abs(np.min(elec_x) - np.max(elec_x))
                # longest dipole calculation available in R2 class
                zlim=[min(self.elec[:,2])-doiEstimate,max(self.elec[:,2])]
        except:
            if xlim==None:
                xlim=[min(self.node[:,0]),max(self.node[:,0])]
            if zlim==None:
                zlim=[min(self.node[:,2]),max(self.node[:,2])]

        if np.diff(xlim) == 0: # protection against thin axis margins 
            xlim=[xlim[0]-2,xlim[1]+2]
        if np.diff(zlim) == 0:
            zlim=[zlim[0]-2,zlim[1]+2]
                
        ##plot mesh! ##
        #compile mesh coordinates into polygon coordinates  
        nodes = np.array([self.node[:,0],self.node[:,2]]).T
        connection = self.connection.copy() # connection matrix 
        if maxDepth is not None: 
            depths = np.array(self.computeElmDepth())
            ikeep = depths < maxDepth
            #truncate connection matrix and plotted array
            connection = connection[ikeep,:]
            X = X[ikeep]
            
        #compile polygons patches into a "patch collection"
        coordinates = nodes[connection]
        if vmin is None:
            vmin = np.min(X)
        if vmax is None:
            vmax = np.max(X)
        
        if edge_color == None or edge_color=='none' or edge_color=='None':
            edge_color='face'#set the edge colours to the colours of the polygon patches

        if contour is False:
            if attr == 'region': # so the default material
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
                xc = self.elmCentre[ikeep,0]
                yc = self.elmCentre[ikeep,2]
            else:
                xc = self.elmCentre[:,0]
                yc = self.elmCentre[:,2]
            zc = np.array(X)
            
            # check for 0 in sigma log
            if attr == 'Sigma_imag(log10)':
                ie = zc != 0
                zc = zc[ie]
                xc = xc[ie]
                yc = yc[ie]
            x = np.array(self.node[:,0])
            y = np.array(self.node[:,2])
            
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
            else: # fallback mode with tricontourf and cropSurface() (topo based on centroids)             
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
            if 'Sensitivity(log10)' not in self.df.keys():
                print('ERROR: No sensitivity attribute found')
            else:
                try:
                    weights = np.array(self.df['Sensitivity(log10)']) #values assigned to alpha channels 
                    if maxDepth is not None:
                        weights = weights[ikeep]
                    if sensPrc is None:
                        thresh = np.log10(0.001*(10**np.nanmax(weights)))
        #                    thresh = np.percentile(weights, 50, interpolation='nearest')
                        x = np.sort(weights)
                        i = np.where(x > thresh)[0][0]
                        x = np.argsort(weights)
                        alphas = np.zeros(self.numel)
                        alphas[:i] = np.linspace(1, 0, len(alphas[:i]))
                        raw_alpha = np.ones((self.numel,4),dtype=float) #raw alpha values 
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
                        x = np.sort(weights)
                        i = np.where(x > thresh)[0][0]
                        x = np.argsort(weights)
                        alphas = np.zeros(self.numel)
                        alphas[:i] = np.linspace(1, 0.2, len(alphas[:i]))
                        raw_alpha = np.ones((self.numel,4),dtype=float) #raw alpha values 
                        raw_alpha[:, -1] = alphas
                        alpha_map = ListedColormap(raw_alpha) # make a alpha color map which can be called by matplotlib
                        #make alpha collection
                        alpha_coll = PolyCollection(coordinates, array=weights, cmap=alpha_map, edgecolors='none', linewidths=0)#'face')
                        #*** the above line can cuase issues "attribute error" no np.array has not attribute get_transform, 
                        #*** i still cant figure out why this is because its the same code used to plot the resistivities 
                        ax.add_collection(alpha_coll)
                    
                except Exception as e:
                    print('Error in the sensitivity overlay:', e)
        
        if electrodes: #try add electrodes to figure if we have them 
            try: 
                x = self.elec.copy()
                if self.iremote is not None: # it's None for quad mesh
                    x = x[~self.iremote, :]
                elecColor = 'ko' if darkMode is False else 'wo'
                ax.plot(x[:,0], x[:,2], elecColor, markersize=4)
            except AttributeError:
                # print("no electrodes in mesh object to plot")
                pass

        # adding interactive display when mouse-over
        centroids = np.array([self.elmCentre[:,0], self.elmCentre[:,2]]).T
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
            
        # decide which attribute to plot, we may decide to have other attritbutes! 
        if attr not in self.df.keys():
            attr = None
        if attr is None:
            attr = self.df.keys()[0]

        X = np.array(self.df[attr])
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
                pvcontour=[],
                pseudo3DContour=False,
                pvdelaunay3d=False,
                pvshow=True,
                darkMode=False,
                cell_picking=False,
                clipping=True):
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
            Axis z limits as `(zmin, zmax)`. 
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
        pvcontour : list of float, optional
            Values of the isosurface to be plotted.
        pseudo3DContour: bool, optional
            If 'True', pseudo 3D plots will be contoured. Only use in case of pseudo 3D surveys.
        pvdelaunay3d : bool, optional
            If `True` a "Delaunay 3D" triangulation filter will be applied on the mesh.
        pvshow : bool, optional
            If `False`, that will prevent calling the `pyvista.Plotter.show()`.
            This is useful in case of subplots.
        darkmode: bool, optional
            Alters coloring of pyvista plot for a darker appearance  
        cell_picking: bool, optional 
            Interactive picking flag, for the use within the UI only. Leave as 
            False. 
        clipping: bool, optional 
            Flag to clip mesh in pyvista window, if tank type problem is expected 
            set to False else set to True. Default is True. Note if False then the 
            X Y and Z limit extents will be ignored. 

        Returns
        -------
        figure : matplotlib figure, pyvista plotter  
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
        # print('+++', self.node[:,1])
        if not isinstance(color_map,str):#check the color map variable is a string
            raise NameError('color_map variable is not a string')
        
        # decide which attribute to plot, we may decide to have other attritbutes! 
        keys = self.df.keys()
        if attr not in keys:
            if attr is not None:
                warnings.warn('Chosen attribute not found in mesh class')
            attr = None
        if attr is None:
            if 'region' in keys:
                attr = 'region'
            else: 
                attr = keys[0]

        X = np.array(self.df[attr])
        color_bar_title = attr
        
        S = np.array([0]*self.numel) 
        if sens:
            try:
                S = np.array(self.sensitivities)
            except AttributeError:
                sens = False
                #print('no sensitivities to plot')
                    
        # determine the depth extent of the survey for bounding box
        if zlim == None:
            zlim = [min(self.node[:,2])-1, max(self.node[:,2])+1]
        
        #work out lateral extents         
        try:
            if self.elec is not None:
                if self.iremote is not None: 
                    iremote = self.iremote
                else:
                    nelec = self.elec.shape[0]
                    iremote = np.array([False]*nelec,dtype=bool)
                    
                if xlim == None:
                    elecx = self.elec[:,0][~iremote]
                    xlim = [min(elecx), max(elecx)]
                if ylim == None:
                    elecy = self.elec[:,1][~iremote]
                    ylim = [min(elecy), max(elecy)]            
        except (AttributeError,TypeError): #if no electrodes present use the node limits 
            pass
        
        if xlim == None: # if still none use the node limits 
            xlim = [min(self.node[:,0])-1, max(self.node[:,0])+1]
        if ylim == None:
            ylim = [min(self.node[:,1])-1, max(self.node[:,1])+1]
          
        # protection against thin axis margins 
        if abs(xlim[0] - xlim[1]) < 0.001:
            xlim=[xlim[0]-2,xlim[1]+2]
        if abs(zlim[0] - zlim[1]) < 0.001:
            zlim=[zlim[0]-2,zlim[1]+2]
        if abs(ylim[0] - ylim[1]) < 0.001:
            ylim=[ylim[0]-2,ylim[1]+2]
        
        # determine min and max of X array 
        if vmin is None:
            vmin = np.min(X)
        if vmax is None:
            vmax = np.max(X)
                    
        #### use pyvista for 3D viewing ####         
        if use_pyvista and pyvista_installed:
            nmesh = self.copy() # this returns a mesh copy
            if clipping:
                nmesh = nmesh.truncateMesh(xlim, ylim, zlim) # truncating is the nly way to reduce the grid extent
            X = nmesh.df[attr].values # NEEDED as truncating change element index
            nmesh.df = pd.DataFrame(X, columns=[color_bar_title]) # make the attr of interest the only attribute            
            
            #save to temperory directory 
            folder = tempfile.TemporaryDirectory()
            fname = os.path.join(folder.name, '__to_pv_mesh.vtk')
            nmesh.vtk(fname)
            self.pvmesh = pv.read(fname) # read in temporary file 
            folder.cleanup()
            
            #### possible to iniate pyvista mesh directly with vertices and connection matrix #### 
            #### however the below code is apparently unstable on pyvista 0.29.0 #### 
            # self.pvmesh = pv.PolyData(self.node, 
            #                           self.connection)
            # self.pvmesh[color_bar_title] = X
                        
            if edge_color is None or edge_color=='none' or edge_color=='None':
                edges = False # then dont show element edges 
            else:
                edges = True
            
            # make a plotter object if not already given or check it
            if ax is None: 
                # ax = BackgroundPlotter()
                ax = pv.Plotter()
                ax.background_color = background_color
            else: # check the ax argument is for pyvista not matplotlib 
                typ_str = str(type(ax))
                if typ_str.find('pyvista') == -1:
                    raise Exception('Error plotting with pyvista, show3D (meshTools.py) expected a pyvista plotter object but got %s instead'%typ_str)
                ax.set_background(background_color)
            
            # Delaunay 3D
            if pvdelaunay3d:
                self.pvmesh = self.pvmesh.cell_data_to_point_data()
                self.pvmesh = self.pvmesh.delaunay_3d()
                color_map = plt.cm.get_cmap(color_map, 14) # subdividing colorbar so it look more like contouring!
            
            # apply threshold
            if pvthreshold is not None:
                if isinstance(pvthreshold, list):
                    if pvthreshold[0] is None:
                        pvthreshold[0] = np.nanmin(X)
                    if pvthreshold[1] is None:
                        pvthreshold[1] = np.nanmax(X)
                self.pvmesh = self.pvmesh.threshold(value=pvthreshold, scalars=attr)
            
            # create isosurfaces
            if len(pvcontour) > 0:
                self.pvmesh = self.pvmesh.cell_data_to_point_data()
                self.pvmesh = self.pvmesh.contour(isosurfaces=pvcontour, scalars=attr)
                
            # set colors for dark mode 
            tcolor = 'k'
            if darkMode:
                tcolor = 'w'
                # elec_color = 'w'
                ax.set_background((0.2,0.2,0.2))
            
             # create pseudo 3D contour plots - only for pseudo 3D surveys
            if pseudo3DContour:
                self.pvmesh = self.pvmesh.cell_data_to_point_data()
                color_map = plt.cm.get_cmap(color_map, 14)

            # show grid
            if pvgrid:
                ax.show_grid(color=tcolor)
            
            # clip mesh to bounding box ... we crop the mesh (only way to reduce its bounds for better viewing - must be after all pvmesh works)
            if clipping:
                self.pvmesh = self.pvmesh.clip_box((xlim[0],xlim[1],ylim[0],ylim[1],zlim[0],zlim[1]),invert=False)
            
            # plot slices or entire mesh
            if np.sum([len(a) for a in pvslices]) > 0: # we have slices
                ax.add_mesh(self.pvmesh.outline(), color=tcolor)
                for i, ss in enumerate(pvslices): # X, Y then Z slice
                    normal = np.zeros(3)
                    normal[i] = 1
                    for s in ss:
                        if ((s > np.nanmin(self.elmCentre[:,i])) &
                            (s < np.nanmax(self.elmCentre[:,i]))):
                            origin = np.zeros(3)
                            origin[i] = s
                            mesh_slice = self.pvmesh.slice(normal=normal, origin=origin)
                            if mesh_slice.number_of_points > 0 or mesh_slice.number_of_cells > 0:
                                ax.add_mesh(mesh_slice,
                                            scalars = attr,
                                            cmap=color_map,
                                            clim=[vmin, vmax],
                                            show_scalar_bar=color_bar,
                                            show_edges=edges,
                                            opacity=alpha,
                                            scalar_bar_args={'color':tcolor})
                            else:
                                print('empty mesh')
            else:        
                if self.pvmesh.number_of_points > 0 or self.pvmesh.number_of_cells > 0:
                    ax.add_mesh(self.pvmesh,
                                scalars = attr, # specify name of vector 
                                cmap=color_map, #matplotlib colormap 
                                clim=[vmin,vmax], #color bar limits 
                                show_scalar_bar=color_bar,#plot the color bar? 
                                show_edges=edges, #show edges
                                # opacity=alpha,
                                scalar_bar_args={'color':tcolor,# 'interactive':True,
                                                 'vertical':False})#,
                                                 #'title_font_size':16,
                                                 #'label_font_size':14})
                else:
                    print('empty mesh')
            
             # then add the electrodes to the plot 
            if electrodes and self.elec is not None: #try add electrodes to figure if we have them 
                try:
                    points = self.elec.copy()[~iremote,:]
                    pvelec = pv.PolyData(points)
                    ax.add_mesh(pvelec, color=elec_color, point_size=10.,
                                render_points_as_spheres=True)
                except AttributeError as e:
                    print("Could not plot 3d electrodes, error = "+str(e))
            
            # show mesh
            if pvshow:
                ax.show()
            # NOTE: because we've truncated the copy of the mesh that pyvista
            # reads in, it should scale correctly in showMesh() even if the
            # mesh object represents the entire mesh with the coarse region.
            
            if cell_picking:
                ax.enable_cell_picking(mesh=self.pvmesh,start=True,through=False)
                #data can then be probed with 
                #data = ax.picked_cells.cell_arrays
                #a = data['name of attribute'] # where a is a numpy array 
                
            
            return ax # exit function
        #### else fall back on to matplotlib scheme ####    
        
        
        #if prism mesh use that function instead 
        if self.ndims==2:
            warnings.warn("Its recommended to use mesh.show() for 2D meshes, results of 3D could be unstable")
        elif self.type2VertsNo() == 6: # use column mesh show instead 
            self.showPrismMesh(color_map = color_map, #displays prism mesh using matplotlib
                color_bar = color_bar,
                xlim = xlim,
                ylim = ylim,
                zlim = zlim, 
                ax = ax,
                electrodes = electrodes,
                sens = sens,
                edge_color = edge_color,
                alpha = alpha,
                vmax = vmax,
                vmin = vmin,
                attr = attr)
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
        
        tmesh = self.truncateMesh(xlim,ylim,zlim)
        
        # search through each element to see if it is on the edge of the mesh, 
        # this step is important as it is very expensive to plot anything in 3D using matplotlib 
        # triangles on the edge of the mesh will be used only once
        
        tmesh.computeNeigh()
        fcon, idx = mc.faces3d(tmesh.connection, 
                               tmesh.neigh_matrix)
        

        node_x = tmesh.node[:,0]
        node_y = tmesh.node[:,1]
        node_z = tmesh.node[:,2]
        
        face_list = [(0,0,0)]*fcon.shape[0]
        for i in range(fcon.shape[0]):
            fnodes = fcon[i,:]
            vert1 = (node_x[fnodes[0]], node_y[fnodes[0]], node_z[fnodes[0]])
            vert2 = (node_x[fnodes[1]], node_y[fnodes[1]], node_z[fnodes[1]])
            vert3 = (node_x[fnodes[2]], node_y[fnodes[2]], node_z[fnodes[2]])
            face_list[i] = (vert1,vert2,vert3)
        
        assign = X[idx]
        sensi = S[idx]
          
        polly = Poly3DCollection(face_list,linewidth=0.5) # make 3D polygon collection
        polly.set_alpha(alpha)#add some transparancy to the elements
        try:
            polly.set_array(np.array(assign))
        except MemoryError:#catch this error and print something more helpful than matplotlibs output
            raise MemoryError("Memory access violation encountered when trying to plot mesh, \n please consider truncating the mesh or display the mesh using paraview.")

        if edge_color != 'face':#attempt to set edge colour 
            try:
                polly.set_edgecolor(edge_color)
            except AttributeError:  #seems to be a bug with matplotlib where this function doesnt work in 3.1
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
            
        if electrodes and self.elec is not None: #try add electrodes to figure if we have them 
            try: 
                ax.scatter(self.elec[:,0],self.elec[:,1],zs=np.array(self.elec[:,2]),
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
        values = np.array(self.df[attr])        
        dimDico = {'x':0,'y':1,'z':2}
        dim = dimDico[axis]
        elms = self.connection.copy()
        nodes = self.node.copy()
        sliceMesh(nodes, elms, values, label=attr, dim=dim, vmin=vmin, vmax=vmax, ax=ax)
        
        
    def show_prism_mesh(self,*args):
        warnings.warn("show_prism_mesh is depreciated, use either show3D or showPrismMesh instead")
        self.showPrismMesh(*args)
        
    def showPrismMesh(self,color_map = 'Spectral',#displays the mesh using matplotlib
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
        """Shows a 3D prism mesh. 
        
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
            
        # decide which attribute to plot, we may decide to have other attritbutes! 
        keys = self.df.keys()
        if attr not in keys:
            attr = None
            warnings.warn('Chosen attribute not found in mesh class or none was chosen')
        if attr is None:
            if 'region' in keys:
                attr = 'region'
            else: 
                attr = keys[0]

        X = np.array(self.df[attr])
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
            zlim=[min(self.node[:,2]), max(self.node[:,2])]
        if xlim==None:
            xlim=[min(self.node[:,0]), max(self.node[:,0])]
        if ylim==None:
            ylim=[min(self.node[:,1]), max(self.node[:,1])]
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
        S = np.array([0]*self.numel) 
        if sens:
            try:
                S = self.sensitivities
            except AttributeError:
                sens = False
                print('no sensitivities to plot')
        
        #call cython code             
        if self.neigh_matrix is None:
            self.computeNeigh()
        neigh = self.neigh_matrix
        face_list, idx = mc.facesPrism(self.connection, self.node, neigh)
        sensi = S[idx]
        assign = X[idx]
            
        #add patches to 3D figure 
        polly = Poly3DCollection(face_list,linewidth=0.5) # make 3D polygon collection
        polly.set_alpha(alpha)#add some transparancy to the elements
        try:
            polly.set_array(np.array(assign))
        except MemoryError:#catch this error and print something more helpful than matplotlibs output
            raise MemoryError("Memory access voilation encountered when trying to plot mesh, \n please consider truncating the mesh or display the mesh using paraview.")
        
        if edge_color != 'face':
            try:
                polly.set_edgecolor(edge_color)
            except AttributeError: # can throw an error with certian versions of matplotlib
                pass
            
        polly.set_cmap(color_map) # set color map 
        polly.set_clim(vmin=vmin, vmax=vmax) # reset the maximum limits of the color map 
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
                ax.scatter(self.elec[:,0],self.elec[:,1],zs=np.array(self.elec[:,2]),
                           s=20, c='k', marker='o')
            except AttributeError as e:
                print("could not plot 3d electrodes, error = "+str(e))
                
        print('Mesh plotted in %6.5f seconds'%(time.time()-t0))
        
        
    #%% 3D selection (designed to be used with the UI only really)
    def _select3Dsurface(self):
        #return surface for 3D selection in forward modelling 
        if self.ndims != 3:
            raise Exception('Only use case for this is function is with 3D meshes')
        if 'gmsh_phys_entity' in self.df.keys():
            a = self.df['gmsh_phys_entity'].values            
            nmesh = self.filterIdx(a==1)
            return nmesh.extractSurface(post_neigh_check=False)
        else:
            try: # cut down to just the electrode area 
                x = self.elec[:,0]
                y = self.elec[:,1]
                tmesh = self.truncateMesh([min(x),max(x)],
                                          [min(y),max(y)])
            except:
                tmesh = self.copy()
            return tmesh.extractSurface()
        
    def pick3Dbox(self,ax=None,xlim=None,ylim=None,zlim=None,electrodes=True, darkMode=False):
        """Pick out regions on the mesh using the box widget (intended for 3D
        forward modelling)

        Parameters
        ----------
        ax : class, optional
            pyvista plotting object. The default is None.
        xlim : tuple, list, array like, optional
            x limit of the selection area. The default is None 
            and assigned automatically.
        ylim : tuple, list, array like, optional
            as with xlim. The default is None.
        zlim : tuple, list, array like, optional
            as with xlim. The default is None.
        electrodes : bool, optional
            If True the electrodes are also plotted on the mesh. The 
            default is True.
        darkmode: bool, optional
            Alters coloring of pyvista plot for a darker appearance 

        Raises
        ------
        Exception
            If mesh is 2D

        Returns
        -------
        pyvista clip box handle
        """
        if self.ndims != 3:
            raise Exception('Only use case for this is function is with 3D meshes')
        if self.type2VertsNo() != 4:
            raise Exception('Function only works on a tetrahedral mesh')
        if ax is None:
            #make pyvista plotter 
            ax = pv.Plotter()
        
        # determine the extent of the survey for bounding box
        if zlim == None:
            zlim = [min(self.node[:,2])-1, max(self.node[:,2])+1]
        
        #work out lateral extents        
        try:
            if self.elec is not None:
                if self.iremote is not None: 
                    iremote = self.iremote
                else:
                    nelec = self.elec.shape[0]
                    iremote = np.array([False]*nelec,dtype=bool)
                elecx = self.elec[:,0][~iremote]
                elecy = self.elec[:,1][~iremote]
                dx = max(elecx) - min(elecx)
                dy = max(elecy) - min(elecy)
                d = np.sqrt(dx**2 + dy**2)
                if xlim == None:
                    xlim = [min(elecx)-d, max(elecx)+d]
                if ylim == None:
                    ylim = [min(elecy)-d, max(elecy)+d]
        except (AttributeError,TypeError): #if no electrodes present use the node limits 
            pass
        if xlim == None:
            xlim = [min(self.node[:,0])-1, max(self.node[:,0])+1]
        if ylim == None:
            ylim = [min(self.node[:,1])-1, max(self.node[:,1])+1]
                
        # if 'elm_id' not in self.df.keys():
        self.df.loc[:,'elm_id'] = np.arange(1,self.numel+1,1)
        
        # surface mesh
        if self.surfaceMesh is None: # do this to avoid recomputing surface mesh
            self.surfaceMesh = self._select3Dsurface() # surface of mesh 
        self.surfaceMesh.show3D(ax=ax,color_map='Greys', color_bar=False,
                     edge_color=None, alpha=0.8, darkMode=darkMode)
        
        # pyvista mesh
        folder = tempfile.TemporaryDirectory()
        fname = os.path.join(folder.name, '__to_pv_mesh.vtk')
        self.vtk(fname)
        pvmesh = pv.read(fname)
        folder.cleanup()
        
        #add all regions != 1 
        reg = self.df['region'].values
        idx = reg>1
        if not all(idx==False):#then add regions > 1. So that background region is transparent
            nmesh = self.filterIdx(idx)
            nmesh.show3D(ax=ax,vmin=1)
            
        #add box clip function
        pvmesh = pvmesh.clip_box((xlim[0],xlim[1],ylim[0],ylim[1],zlim[0],zlim[1]),invert=False)
        ax.add_mesh_clip_box(pvmesh, color='white',opacity=0.4)
        
        if self.elec is not None and electrodes: #try add electrodes to figure if we have them 
            points = self.elec.copy()
            elec_color = 'k' # give option to change this in future? 
            pvelec = pv.PolyData(points)
            ax.add_mesh(pvelec, color=elec_color, point_size=10.,
                        render_points_as_spheres=True)
        
        return ax.box_clipped_meshes[0]
    
    def addRegion3D(self, clipped):
        if self.ndims != 3:
            raise Exception('Only use case for this function is with 3D meshes')
        idx = np.unique(np.asarray(clipped.cell_arrays['elm_id'], dtype=int))-1
        self.df.loc[idx, 'region'] = np.max(self.df['region']) + 1 
        
        in_elem = np.array([False]*self.numel,dtype=bool)
        in_elem[idx] = True
                
        nmesh = self.filterIdx(in_elem)
        
        return nmesh
            
        
    #%% mesh zoning     
    def assignZone(self,poly_data):
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
            
        numel=self.numel#number of elements 
        elm_xz=np.array([self.elmCentre[:,0],self.elmCentre[:,2]]).T#centriods of 2d mesh elements 
        material_no=np.zeros(numel,dtype=int)#attribute number
        
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
        self.addAttribute(material_no,'zone')
        return material_no     
        
    def assignZoneAttribute(self,attr_list,new_key,zone=None):
        """ Asssigns values to the mesh which depend on region / material only. E.G 
        a single resistivity value.
            
        Parameters
        ----------
        attr_list : list
            Values corresponding to a material number in the mesh. eg. if you had 3 regions in the mesh then you give
            `[resistivity1,resistivity2,resistivity3]`.
        new_key : string
            Key identifier assigned to the attribute in the df. 
        zone : array or list
            Integers starting at 0 or 1, and ascend in intervals of 1, which 
            correspond to a material in the mesh returned from assign_attr_ID.
            Should have the same length as the number of elements in the mesh.
        
        Notes  
        -----
        Mesh object will now have the new attribute added once the function is run.
        Use the `mesh.show()` (or `.draw()`) function to see the result. 
        """ 
        if zone is None:
            zone = self.zone
            
        if len(zone) != self.numel:
            raise ValueError("Mismatch between the number of elements and material propeties")
        
        new_para=np.array([0]*self.numel)
        
        if min(zone)==1:#cor_fac allows for compatability with an index system starting at 1 or 0 
            cor_fac=1
        else:
            cor_fac=0
            
        zone = np.array(zone) - cor_fac
        
        for i in range(len(attr_list)):
            idx = zone == i
            new_para[idx] = attr_list[i]
        
        self.df[new_key] = new_para
        self.no_attributes += 1            

    def applyFunc(self,mesh_paras,material_no,new_key,function,*args):
        """ Applies a function to a mesh by zone number and mesh parameter.
        
        Parameters
        ----------
        mesh_paras : array like
            Mesh parameters from which new parameters are calculated.
        material_no : array like of ints
            Material type assigned to each element, should be numbered consectively from 1 to n. in the form 1 : 1 : 2 : n.
            ...ie if you have 2 materials in the mesh then pass an array of ones and twos. zeros will be ignored. 
        new_key : string
            Key assigned to the parameter in the df. DOES NOT default.
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
        new_para=[0]*self.numel
        #iterate through each set of argument variables
        for iteration in range(len(args[0])):
            parameters=[items[iteration] for items in args]#return parameters 
            parameters.insert(0,0)#this adds an element to the front of the parameters which can be swapped out to mesh_paras
            for i in range(self.numel):
                if material_no[i]==iteration+1:#does the material match the iteration? 
                    parameters[0]=mesh_paras[i]#change parameter value at start of variables list
                    new_para[i]=function(*parameters)#compute new parameter   
        self.df[new_key] = new_para
        self.no_attributes += 1
        #return new_para
        
    def computeElmDepth(self):
        """ Compute the depth of elements relative to the surface of the mesh.
        """
        if self.ndims == 2: # use 1D interpolation
            xz = self.extractSurface()
            datum_x = xz[0]
            datum_z = xz[1]
            
            elm_x = np.array(self.elmCentre[:,0])
            elm_z = np.array(self.elmCentre[:,2])
            min_idx = np.argmin(datum_x)
            max_idx = np.argmax(datum_x)
            Z = np.interp(elm_x,datum_x,datum_z,left=datum_z[min_idx],
                          right=datum_z[max_idx])
            depth = Z - elm_z
            self.df['depths'] = depth
            self.no_attributes += 1
            return depth
        elif self.ndims == 3: # use 2D interpolation
            elm_x = np.array(self.elmCentre[:,0])
            elm_y = np.array(self.elmCentre[:,1])
            elm_z = np.array(self.elmCentre[:,2])  
            #use interpolation to work out depth to datum 
            surmesh = self.extractSurface()
            Z = interp.triangulate(elm_x, elm_y, surmesh.node[:,0], 
                                   surmesh.node[:,1], surmesh.node[:,2])
            depth = Z - elm_z
            self.df['depths'] = depth # add cell depths to attribute cache
            self.no_attributes += 1
            return depth
        
    #%% Gradients 
    def downslopeID(self,attr='Resistivity'):
        """
        Get the index of 'downslope' elements for a given attribute. 

        Parameters
        ----------
        attr : string, optional
            Key inside of mesh dataframe. The default is 'Resistivity'.

        Returns
        -------
        idx : nd array 
            Indices of downslope elements. 
        """
        
        #get neighbour matrix  
        if self.neigh_matrix is None:
            self.computeElmDepth()
        neigh = self.neigh_matrix
        
        Z = self.df[attr].values
        
        #find the steepest gradient for each neighbour 
        Zcomp = np.zeros(3)
        idx = np.zeros(self.numel,dtype=int)
        for i in range(neigh.shape[0]):
            for j in range(neigh.shape[1]):
                ni = neigh[i,j]
                if ni != -1:
                    Zcomp[j] = Z[i] - Z[ni]
            #where most postive the Z value is downslope
            arg = np.argmax(Zcomp)
            idx[i] = neigh[i,arg]
            
            #filter out instances where downslopeID is actually upslope  
            check = Z[i] - Z[idx[i]]
            if check < 0:
                idx[i] = -1 # minus means no downslope cell (like edge cases )
        
        return idx 
    
    def node2ElemAttr(self, node_array, name = 'NodeArray'):
        """
        Maps node attributes onto mesh elements 

        Parameters
        ----------
        node_array : array like 
            Array with same number of entries as there are nodes in the mesh. 
        name : str, optional
            Name of node attribute, will be used to index array in mesh dataframe. 
            The default is 'NodeArray'.

        Returns
        -------
        arr : nd array 
            Node array mapped on mesh elements .

        """
        if len(node_array) != self.numnp:
            raise ValueError("Length of node array does not match the number of nodes inside the mesh")
            
        node_arr = np.array(node_array) # array of node values 
        arr = np.mean(node_arr[self.connection], axis = 1) # values of array on the node values 
        self.addAttribute(arr,name)
        return arr 
    
    #%% electrode movements         
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
        
        has_nodes = True
        if self.eNodes is None:
             has_nodes = False
    
        if has_nodes and len(new_z) != len(self.eNodes):
            raise ValueError("mismatch between the length of new electrode position array and old array")
            
        if new_y is None:
            new_y = np.zeros_like(new_x)
            
        node_x = self.node[:,0]
        node_y = self.node[:,1]
        node_z = self.node[:,2]
            
        node_in_mesh = [0]*len(new_x)    
        for i in range(len(new_x)):
            sq_dist = (node_x - new_x[i])**2 + (node_y - new_y[i])**2 + (node_z - new_z[i])**2 # find the minimum square distance
            node_in_mesh[i] = np.argmin(sq_dist) # min distance should be zero, ie. the node index.
            if has_nodes and debug:
                if node_in_mesh[i] != self.eNodes[i]:
                    print("Electrode %i moved from node %i to node %i"%(i,node_in_mesh[i],self.eNodes [i]))#print to show something happening
        
        self.setElecNode(node_in_mesh) # update e_node parameter
        if len(np.unique(node_in_mesh)) != len(new_x):
            warnings.warn("The number of new electrode nodes does not match the number of electrodes, which means a duplicated node is present! Please make mesh finer.")   
     
        return np.array(node_in_mesh, dtype=int) # note this is the node position with indexing starting at 0. 

    #%% write mesh to file 
    def findIdirichlet(self):
        """Find the best node for the dirchlet node 
        Returns
        -------
        idirchlet: int 
            A node far away as possible from each of the electrodes on the boundary of the 
            mesh. (Add one if using inside of mesh(3d).dat)
        """
        if self.eNodes is None:
            #this function cannot run because it assumes that it knows the electrode nodes
            warnings.warn('Taken last node of mesh as dirchlet node, as no electrode nodes assigned to mesh class')
            return self.numnp-1  #in which case take last node as idirchlet 
        
        if self.neigh_matrix is None:
            self.computeNeigh()
        neigh = self.neigh_matrix.copy()
        edge = self.connection[np.min(neigh,axis=1) == -1,:]
        edge_nodes = self.node[np.unique(edge.flatten()),:]
        elec = self.node[self.eNodes,:]
        # find average electrode position 
        e = np.c_[np.mean(elec[:,0]),np.mean(elec[:,1]),np.mean(elec[:,2])]
        
        idirichlet = np.argmax(np.sqrt(np.sum((edge_nodes - e)**2, axis=1)))
        
        return idirichlet

        
    def write_dat(self, file_path='mesh.dat'):
        warnings.warn('write_dat is depreciated, use dat instead')
        self.dat(file_path)
        
    def dat(self, file_path='mesh.dat'):
        """Write a mesh.dat kind of file for mesh input for R2/R3t. 
        
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
            idirichlet = self.findIdirichlet() + 1 
            if self.ndims == 3:
                fid.write('%i %i %i 0 %i\n'%(self.numel,self.numnp,1,self.type2VertsNo()))
            else:
                fid.write('%i %i %i\n'%(self.numel,self.numnp,idirichlet))
            
            try:
                param = np.array(self.df['param'])
            except:
                param = np.arange(1,self.numel+1)
            try:
                zone = np.array(self.df['zones'])
            except:
                zone = np.ones(self.numel,dtype=int)

            #write out elements         
            no_verts = self.type2VertsNo()
            for i in range(self.numel):
                elm_no=i+1
                fid.write("%i "%elm_no)
                [fid.write("%i "%(self.connection[i][k]+1)) for k in range(no_verts)]
                fid.write("%i %i\n"%(param[i],zone[i]))
    
            #now add nodes
            if self.ndims == 3:
                nidx = [0,1,2]
            else:
                nidx = [0,2]
                
            for i in range(self.numnp):
                fid.write('{:<16d} '.format(i+1)) # node number 
                [fid.write('{:<16.8f} '.format(self.node[i,k])) for k in nidx] # node coordinates 
                fid.write('\n')#drop down a line 
                    
    def datAdv(self, file_path='mesh.dat', iadvanced=True):
        """Write a mesh.dat kind of file for mesh input for R2/R3t. Advanced format
        which includes the neighbourhood and conductance matrix. 
        
        Parameters
        ----------
        file_path : str, optional
            Path to the file. By default 'mesh.dat' is saved in the working directory.
        iadvanced : bool, optional
            If `True`, use the advanced mesh format if `False`, use the normal
            mesh format.
        """
        if not isinstance(file_path,str):
            raise TypeError("Expected string argument for file_path")
        if self.ndims == 2:
            raise TypeError('Advanced mesh format not avialable with 2D meshes currently')
            
        # find furthest node from first electrode to be dirichlet node
        idirichlet = self.findIdirichlet() + 1 
            
        # element parameters 
        param = np.asarray(self.df['param'],dtype=int)
        try: 
            zone = np.asarray(self.df['zones'],dtype=int)
        except:
            zone = np.ones(self.numel,dtype=int)
            
            
        adv_flag = int(iadvanced)
        #compute neighbourhood matrix 
        if self.neigh_matrix is None:
            self.computeNeigh()
        neigh = self.neigh_matrix.copy()
        neigh = mc.sortNeigh(neigh) # organise for (c)R3t input  
        neigh += 1 # add one for FOTRAN indexing 
        
        if self.NsizeA is None:#then the finite element conductance matrix needs calculating 
            self.computeNconnec()
            
        # NsizeA = self.NsizeA
        fconm = self.fconm.copy() + 1 #(add 1 for FORTRAN indexing)
        ### write data to mesh.dat kind of file ###
        # open mesh.dat for input      
        with open(file_path, 'w') as fid:
            # write to mesh.dat total num of elements and nodes
            if self.ndims == 3:
                fid.write('%i\t%i\t%i\t%i\t%i\t%i\n'%(self.numel,self.numnp,1,0,self.type2VertsNo(),adv_flag)) # flags 
            else:
                fid.write('%i\t%i\t%i\t%i\n'%(self.numel,self.numnp,idirichlet,adv_flag))
            # write out elements         
            no_verts = self.type2VertsNo()
            for i in range(self.numel):
                fid.write("{:<16d} ".format(i+1)) # add element number 
                [fid.write("{:<16d} ".format(self.connection[i,k]+1)) for k in range(no_verts)] # connection matrix 
                fid.write("{:<16d} {:<16d} ".format(param[i],zone[i])) # parameter and zones 
                # add neighbours in advanced mode
                [fid.write('{:<16d} '.format(neigh[i,k])) for k in range(neigh.shape[1])]#add neighbours 
                fid.write('\n')#drop down a line 
    
            # now add nodes
            if self.ndims == 3:
                nidx = [0,1,2]
            else:
                nidx = [0,2]
                
            for i in range(self.numnp):
                fid.write('{:<16d} '.format(i+1)) # node number 
                [fid.write('{:<16.8f} '.format(self.node[i,k])) for k in nidx] # node coordinates 
                #add conductance matrix in advanced mode 
                [fid.write('{:<16d} '.format(fconm[i,k])) for k in range(fconm.shape[1])] 
                fid.write('\n')#drop down a line 

            if self.ndims == 3: 
                fid.write('{:d}'.format(idirichlet))
                    

    def write_vtk(self, file_path="mesh.vtk", title=None, replace_nan=-9999):
        warnings.warn('write_vtk is depreciated, use vtk instead')
        self.vtk(file_path, title, replace_nan)
        
    def vtk(self, file_path="mesh.vtk", title=None, replace_nan=-9999):
        """Writes a vtk file for the mesh object, everything in the df
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
        fh.write("POINTS %i double\n"%self.numnp)
        node_x = self.node[:,0]
        node_y = self.node[:,1]
        node_z = self.node[:,2]
        
        for i in range(self.numnp):
            fh.write("%16.8f\t%16.8f\t%16.8f\n"%(node_x[i],node_y[i],node_z[i]))
        
        #define the connection matrix    
        no_verts = self.type2VertsNo()
        no_readable = self.numel*(1+no_verts)
        fh.write("CELLS %i %i\n"%(self.numel,no_readable))
        for i in range(self.numel):
            fh.write("%i\t"%no_verts)
            for k in range(no_verts):
                fh.write("{}    ".format(self.connection[i][k]))
            fh.write("\n")
        
        #cell types
        fh.write("CELL_TYPES %i\n"%self.numel)
        [fh.write("%i "%self.cell_type[0]) for i in range(self.numel)];fh.write("\n")
        
        #write out the data
        fh.write("CELL_DATA %i\n"%self.numel)
        df = self.df
        for i,key in enumerate(df.keys()):
            fh.write("SCALARS %s double 1\n"%key.replace(' ','_'))
            fh.write("LOOKUP_TABLE default\n")
            X = np.array(df[key])
            X[np.isnan(X)]=replace_nan
            [fh.write("%16.8f "%X[j]) for j in range(self.numel)]
            fh.write("\n")
        
        #finish writing
        fh.write("POINT_DATA %i\n"%self.numnp)     
        if self.ptdf is not None and len(self.ptdf.keys())>0:
            for i,key in enumerate(self.ptdf.keys()):
                fh.write("SCALARS %s double 1\n"%key.replace(' ','_'))
                fh.write("LOOKUP_TABLE default\n")
                X = np.array(self.ptdf[key])
                X[np.isnan(X)]=replace_nan
                [fh.write("%16.8f "%X[j]) for j in range(self.numnp)]
                fh.write("\n")
        fh.close()
    
    def write_attr(self,attr_key=None,file_name='_res.dat'):
        warnings.warn('write_attr is depreciated use writeAttr instead')
        self.writeAttr(attr_key,file_name)
    
    def writeAttr(self,attr_key=None,file_name='_res.dat'):
        """ Writes a attribute to a _res.dat type file. file_name entered
        seperately because it will be needed for the R2 config file.
        The reason for this function is so you can write a forward model 
        parameter input file. 
        
        Parameters
        ----------
        attr_key: string
            Key identifying the attr to be written in the mesh object df.
        file_name: string, optional
            Name of the _res.dat type file.
        """
        #formality checks
        #double check the fname is not longer than 15 characters as the code doesnt like it apparently 
        if len(ntpath.basename(file_name))>15:
            raise NameError("File name for _res.dat type file cannot be longer than 15 characters")
            
        if attr_key is None: 
            attr_key = 'Resistivity'
        
        if attr_key not in self.df.keys():
            raise NameError('Given key is not inside of the mesh.df')
            
        if isinstance(file_name,str)==False:
            raise NameError("file_name argument must be a string")
        
        #the format of the _res.dat file is such that
        #| x coordinate | y coordinate | value | log(value) | 
        fh = open(file_name,'w')#open file handle 
        x_coords=self.elmCentre[:,0]#get element coordinates
        y_coords=self.elmCentre[:,1]
        z_coords=self.elmCentre[:,2]
        values=self.df[attr_key]
        log_values=np.log10(np.array(values))
        if self.ndims==3:
            for i in range(self.numel):
                fh.write("\t{: 10.5e}\t{: 10.5e}\t{: 10.5e}\t{: 10.5e}\t{: 10.5e}\n".format(x_coords[i],y_coords[i],z_coords[i],values[i],log_values[i]))
        else:
            for i in range(self.numel):
                fh.write("\t{: 10.5e}\t{: 10.5e}\t{: 10.5e}\t{: 10.5e}\n".format(x_coords[i],z_coords[i],values[i],log_values[i]))
            
        fh.close()
        
    def toCSV(self,file_name='mesh.csv'):
        """ Write a .csv file of the mesh, the first 3 columns are the element 
        centres at coordinates x,y,z and the rest of the columns are the 
        attributes in the df

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
        x_coords=self.elmCentre[:,0]#get element coordinates
        y_coords=self.elmCentre[:,1]
        z_coords=self.elmCentre[:,2]
        df = self.df.copy().reset_index()
        keys0= self.df.keys()
        keys=[]
        #ignore the x y z columns if already in df
        ignore = ['X','Y','Z']
        for k in keys0:
            if k not in ignore:
                keys.append(k)
        
        #open file 
        fh = open(file_name,'w')
        fh.write('x,y,z') # write xyz headers 
        [fh.write(','+key) for key in keys] # write attribute headers 
        fh.write('\n') # drop a line 
        for i in range(self.numel):
            line = '{:f},{:f},{:f}'.format(x_coords[i],y_coords[i],z_coords[i])
            for key in keys:
                line += ',{:}'.format(df[key][i])
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
                
        self.vtk(fname)#write vtk to working directory with all associated attributes
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
            
    
    def exportTetgenMesh(self,prefix='mesh',zone=None, debug=False):
        """Export a mesh like the tetgen format for input into E4D. 
        This format is composed of several files. Currently only tested for 
        3D surface array like surveys. Assumes the sides of the mesh are parrallel 
        to the x and y plane and that the surface has gentle to no topography. 
        
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
            if len(zone) != self.numel:
                raise ValueError('Number of zone array elements given to exportTetgenMesh does not match the number of elements in the mesh.')
        else:
            zone = [1]*self.numel # all elements are inside zone 1 
        
        if debug: # print outputs? 
            def stream(s,**kwargs):
                print(s,**kwargs)
        else:
            def stream(s,**kwargs):
                pass 
            
        stream('## Exporting to tetgen format  ##')
            
        stream('Computing boundary conditions ...', end ='')# this an intense process:        
        #NB: faces on side/bottom of mesh are given a marker of 2. Face on top of mesh denoted 1. 
        #NB: Nodes on side of mesh are given a marker of 2, 1 if on top of the mesh and 0 if inside the mesh 

        if self.neigh_matrix is None: # compute neighbour matrix 
            self.computeNeigh() # this will find the element neighbours and the elements which lie on the outside of the mesh! 
            
        ## offload boundary condition calls to cython 
        neigh = self.neigh_matrix
        faces, idx = mc.faces3d(self.connection, neigh) # get outer faces 
        tconnec = self.connection[idx,:] # truncated connection matrix 
        node_bd, face_bd = mc.boundcall(tconnec, faces, self.node) # face and node boundary conditions 
        self.addPtAttribute(node_bd,'node_bd')
        stream('done')
                
        ## output .node file
        stream('writing .node file... ',end='')
        fh = open(prefix+'.1.node','w')
        #header line : 
        #<# of points> <dimension (3)> <# of attributes> <boundary markers (0 or 1)>
        fh.write('{:d}\t{:d}\t{:d}\t{:d}\n'.format(self.numnp,self.ndims,1,1))
        #all other lines: 
        #<point #> <x> <y> <z> [attributes] [boundary marker]
        for i in range(self.numnp):
            line = '{:d}\t{:f}\t{:f}\t{:f}\t{:d}\t{:d}\n'.format((i+1),
                                                        self.node[i,0],
                                                        self.node[i,1],
                                                        self.node[i,2],
                                                        1,
                                                        node_bd[i])
            fh.write(line)
        fh.write('# exported from meshTools module in ResIPy electrical resistivity processing package')
        fh.close()   
        stream('done.')
            
        ## output .ele file 
        stream('writing .ele file... ',end='')
        fh = open(prefix+'.1.ele','w')
        #First line: <# of tetrahedra> <nodes per tet. (4 or 10)> <region attribute (0 or 1)>
        fh.write('{:d}\t{:d}\t{:d}\n'.format(self.numel,self.type2VertsNo(),1))
        #Remaining lines list # of tetrahedra:<tetrahedron #> <node> <node> ... <node> [attribute]
        for i in range(self.numel):
            line = '{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\n'.format((i+1),
                                                         self.connection[i,0]+1,#need to add one because of fortran indexing 
                                                         self.connection[i,1]+1,
                                                         self.connection[i,2]+1,
                                                         self.connection[i,3]+1,
                                                         zone[i])
            fh.write(line)
        fh.write('# exported from meshTools module in ResIPy electrical resistivity processing package')
        fh.close()
        stream('done.')
        
        #output .trn file 
        stream('writing .trn file... ',end='')
        fh = open(prefix+'.trn','w')
        fh.write('0\t0\t0')
        fh.close()
        stream('done.')
        
        #write .face file - which describes elements on the outer edges of the mesh
        
        stream('Writing .face file... ',end='') 

        fh = open(prefix+'.1.face','w')
        #First line: <# of faces> <boundary marker (0 or 1)>
        fh.write('{:d}\t{:d}\n'.format(faces.shape[0],1))#header line 
        for i in range(faces.shape[0]):
            line = '{:d}\t{:d}\t{:d}\t{:d}\t{:d}\n'.format((i+1),
                                                        faces[i,0]+1,
                                                        faces[i,1]+1,
                                                        faces[i,2]+1,
                                                        face_bd[i])
            fh.write(line)
            
        
        fh.write('# exported from meshTools module in ResIPy electrical resistivity processing package')    
        fh.close()
        stream('done.')
        
        ## out .neigh file
        stream('writing .neigh file... ',end='')
        neigh_matrix = neigh.copy()
        neigh_matrix+=1 # add one for tetgen indexing
        neigh_matrix[neigh_matrix==0]=-1#-1 for face elements 
        
       
        fh = open(prefix+'.1.neigh','w') # write to file 
        fh.write('{:d}\t{:d}\n'.format(self.numel,self.type2VertsNo())) # header line         
        elm_id = np.arange(self.numel) + 1
        for i in range(self.numel):
            line = '{:d}\t{:d}\t{:d}\t{:d}\t{:d}\n'.format(elm_id[i],
                                                            neigh_matrix[i,0],
                                                            neigh_matrix[i,1],
                                                            neigh_matrix[i,2],
                                                            neigh_matrix[i,3])
            fh.write(line)
            
        fh.write('# exported from meshTools module in ResIPy electrical resistivity processing package')    
        fh.close()
        stream('done.')
        
    def saveMesh(self, fname, ftype=None):
        """Save mesh into a file. Avaialble formats are .dat, .vtk and .node

        Parameters
        ----------
        fname : TYPE
            DESCRIPTION.
        """
        if not isinstance(fname,str):
            raise TypeError('fname needs to be a string!')
        
        #determine file type 
        atypes = ['dat','node','vtk','csv']
        if ftype is None:#guess
            for a in atypes:
                if fname.endswith('.' + a):
                    ftype = a
                    break
        
        if ftype not in atypes:
            raise NameError('Unregocnised mesh file format, avaialable types are %s'%str(atypes))
        
        #add extension if not already there 
        if not fname.endswith('.'+ftype):
            fname = fname + '.' + ftype
            
        #call save file function 
        if ftype == 'dat':
            self.dat(fname)
        elif ftype == 'vtk':
            self.vtk(fname)
        elif ftype == 'csv':
            self.toCSV(fname)
        elif ftype == 'node':
            self.exportTetgenMesh(fname.replace('.node',''))
        
        
    def writeRindex(self,fname):
        """Write out the neighbourhood matrix for R3t. 

        Parameters
        ----------
        fname : str
            Path to written file.
        """
        if self.type2VertsNo() == 6 or self.type2VertsNo()==8:#not a tetrahedra 3d mesh 
            raise Exception("Sorry neighbour calculation not available yet with this mesh type")
        elif self.type2VertsNo() == 4 and self.ndims==2:#then its a quad mesh 
            raise Exception("Sorry neighbour calculation not available yet with this mesh type")
         
        dim = self.connection.shape[1]

        if self.neigh_matrix is None: 
            self.computeNeigh()
        
        neigh = self.neigh_matrix + 1
        fh = open(fname,'w')
        fh.write('%i\t%i\n'%(self.numel,dim+2))
        for i in range(self.numel):
            line = '{:d}'.format(i+1)
            for j in range(dim):
                line += '\t{:d}'.format(neigh[i,j])
            line += '\t{:d}\n'.format(0)
            fh.write(line)
        fh.close()
        
#%%  Moving mesh
def moveMesh2D(meshObject, elecLocal, elecGrid):
    """Move mesh object to a certain place on the grid.
        X, Y only. No Z relocation.
    
    Parameters
    ----------
    meshObject : class object
        Mesh class object.
    elecLocal: dataframe
        dataframe of electrodes on local 2D grid (i.e., y = 0)
        'remote' column (bool) required.
    elecGrid: dataframe
        dataframe of electrodes on 3D grid (i.e., y != 0)
        'remote' column (bool) required
        
    Returns
    -------
    mesh : class object
        Copy of moved Mesh class object.
    """
    mesh = deepcopy(meshObject) # Mesh.copy() not working here! not sure why!
    xLocal = elecLocal['x'][~elecLocal['remote']].copy().values
    xGrid = elecGrid['x'][~elecGrid['remote']].copy().values
    yGrid = elecGrid['y'][~elecGrid['remote']].copy().values
    if np.all(xGrid == xGrid[0]): # line is vertical - x is constant
        mesh.elec[:,0][~mesh.iremote] = np.ones_like(xLocal) * xGrid[0]
        mesh.elec[:,1][~mesh.iremote] = xLocal + yGrid[np.argmin(yGrid)]
        mesh.node[:,1] = mesh.node[:,0] + yGrid[np.argmin(yGrid)]
        mesh.node[:,0] = np.zeros_like(mesh.node[:,0]) + xGrid[0]

    elif np.all(yGrid == yGrid[0]): # line is horizontal - y is constant
        mesh.elec[:,0][~mesh.iremote] = xLocal + xGrid[np.argmin(xGrid)]
        mesh.elec[:,1][~mesh.iremote] = np.ones_like(xLocal) * yGrid[0]
        mesh.node[:,0] = mesh.node[:,0] + xGrid[np.argmin(xGrid)]
        mesh.node[:,1] = np.zeros_like(mesh.node[:,1]) + yGrid[0]
        
    else: 
        mx = np.polyfit(xLocal,xGrid,1)
        my = np.polyfit(xGrid, yGrid,1)
        mesh.elec[:,0][~mesh.iremote] = mx[0] * xLocal + mx[1]
        mesh.elec[:,1][~mesh.iremote] = my[0] * xGrid + my[1]
        mesh.node[:,0] = mx[0] * mesh.node[:,0] + mx[1]
        mesh.node[:,1] = my[0] * mesh.node[:,0] + my[1]

    return mesh   

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
    
    if vtk_ver.find('vtk')==-1:#version handling 
        raise ImportError("Unexpected file type... ")
    elif vtk_ver.find('3.0')==-1:#not the development version for this code
        warnings.warn("Warning: vtk manipulation code was developed for vtk datafile version 3.0, unexpected behaviour may occur in resulting mesh")
        if vtk_ver.find('4') > -1:
            try:
                # print('Looks like the file is legacy format 4, attempting to import data with newer parser...',end='')
                mesh = vtk_import_fmt4(file_path,order_nodes)
                # print('Success')
                return mesh
            except:
                # print('Failed')
                pass
    
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
        numel=int(elm_info[1])
    except IndexError: # quick bug fix
        elm_info=fid.readline().strip().split()#read line with cell data
        numel=int(elm_info[1])
    
    if numel ==0: 
        raise ImportError("No elements in vtk file to import!")
        
    #read in first element and decide what mesh type it is 
    elm_data=fid.readline().strip().split()
    npere = int(elm_data[0]) # number of vertices per element 
    elm_num = [0]*numel
    con_mat = [[0]*numel for i in range(npere)]
    for j in range(npere):
        con_mat[j][0] = int(elm_data[j+1])
    
    #read in the rest of the elements 
    for i in range(1,numel):
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
            while len(cell_type) < numel:
                cell_type += [float(x) for x in cell_attr_dump[line_idx].split()]
                line_idx += 1
            break
    
    fid.close()
    # read through cell attributes to find the relevant parameter table?
    
    #find scalar values in the vtk file
    num_attr = 0
    attr_dict = {}
    pt_dict = {}
    cells = True 
    lookup = numel # number of possible lookup entries (depends whether looking at cell or point data)

    for i,line_info in enumerate(cell_attr_dump):
        if line_info.find("CELL_DATA") == 0: 
            cells = True
            lookup = numel
        if line_info.find("POINT_DATA") == 0:
            cells = False 
            lookup = numnp 
        if line_info.find("SCALARS") == 0:
            attr_title = line_info.split()[1]
            #check look up table
            if cell_attr_dump[i+1].split()[1] != "default":
                warnings.warn("unrecognised lookup table type")
            line_idx = i+2
            #values=[float(k) for k in cell_attr_dump[i+2].split()]
            array = []
            #while loop to get all scalar values 
            while len(array) < lookup:
                array += [float(x) for x in cell_attr_dump[line_idx].split()]
                line_idx += 1
            if cells:
                attr_dict[attr_title] = array
                num_attr += 1
            else:
                pt_dict[attr_title] = array 
            
    #if cell_type[0] == 5 or cell_type[0] == 8 or cell_type[0] == 9: # then its a 2D mesh
    if title == 'Output from cR2' or title == 'Output from R2': # account for the fact the y and z columns should be swapped 
        mesh = Mesh(node_x,#x coordinates of nodes 
                    node_z,#y coordinates of nodes
                    node_y,#z coordinates of nodes  
                    np.array(con_mat).T,#nodes of element vertices
                    cell_type,#according to vtk format
                    file_path,
                    order_nodes) #nb: nodes should ordered already, not if quad  
    else:
        mesh = Mesh(node_x,#x coordinates of nodes 
                    node_y,#y coordinates of nodes
                    node_z,#z coordinates of nodes 
                    np.array(con_mat).T,#nodes of element vertices
                    cell_type,#according to vtk format
                    file_path,
                    order_nodes) 
    
    #add attributes / cell parameters 
    for key in attr_dict.keys():
        mesh.df[key] = attr_dict[key]
        #see if sensivity output from R2/R3t is inside the .vtk 
        if key in ['Sensitivity_map(log10)','Sensitivity(log10)']:
            mesh.addSensitivity(np.array(attr_dict[key]))
    
    mesh.mesh_title = title
    return mesh

def vtk_import_fmt4(file_path,order_nodes=True):
    """Vtk importer for newer format (not fully validated)
    """
    if os.path.getsize(file_path)==0: # So that people dont ask me why you cant read in an empty file, throw up this error. 
        raise ImportError("Provided mesh file is empty! Check that (c)R2/3t code has run correctly!")
    #open the selected file for reading
    fid = open(file_path,'r')
    
    dump = fid.readlines()
    title = dump[1].strip()
    fid.close()
    
    d = {} # cache of element attributes 
    sid = {}
    #read in header info and perform checks to make sure things are as expected
    for i,line in enumerate(dump):
        if "BINARY" in line:
            raise ImportError("expected ASCII type file format, not binary")
        if "POINTS" in line:
            l = line.split()
            numnp = int(l[1])        
            node_start = i
        if "POLYGONS" in line:
            l = line.split()
            numel = int(l[1])        
            if numel == 0:
                raise ImportError("No elements in vtk file to import!")
            elem_start = i
        if 'FIELD' in line:
            fIEld_id = i 
            
    for i,line in enumerate(dump):
        if "double" in line and i >  fIEld_id:
            l = line.split()
            d[l[0]] = []
            sid[l[0]] = i+1
            
    #now read in node data
    node = np.zeros((numnp,3),dtype=float)
    c = 0
    li = node_start+1
    while c < numnp:
        line = [float(x) for x in dump[li].split()]
        j=0
        i=0
        while j < len(line):
            node[c,i] = line[j]
            i+=1
            if i == 3:
                i=0#reset
                c+=1
            j+=1
        li+=1
        
    #now read in element data
    li = elem_start+1
    line = [int(x) for x in dump[li].split()]
    npere = line[0]
    kx = np.zeros((numel,npere),dtype=float)
    
    for i in range(numel):
        line = [int(x) for x in dump[li].split()]
        for j in range(npere):
            kx[i,j] = line[j+1]
        li+=1
        
    #cell attributes
    for i,key in enumerate(d.keys()):
        li = sid[key]
        while len(d[key]) < numel:
            line = [float(x) for x in dump[li].split()]
            d[key] += line
            c += len(line)
            li+=1
            
    if all(node[:,1]==0) and npere==4:#quad mesh 
        cell_type = [9]*numel
    elif npere ==4:
        cell_type = [10]*numel
    elif npere == 3:
        cell_type = [5]*numel
    elif npere == 6:
        cell_type = [13]*numel
    elif npere == 8:
        cell_type = [11]*numel 
    else:
        raise ImportError('Unsupported cell type in vtk file')
                
    mesh = Mesh(node[:,0],#x coordinates of nodes 
                node[:,1],#y coordinates of nodes
                node[:,2],#z coordinates of nodes 
                kx,#nodes of element vertices
                cell_type,#according to vtk format
                file_path,
                order_nodes) 
    
    #add attributes / cell parameters 
    for key in d.keys():
        mesh.df[key] = np.array(d[key])
    
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
    if len(header)==3: # its a 2D mesh 
        flag_3d = False
        npere = 3
    elif len(header)==5:#old mesh format (e.g. R3t 1.8)
        flag_3d = True
        npere = int(header[-1])
    else:
        flag_3d = True#newer mesh format (R3t v2.20)
        npere = int(header[-2])
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
                node_data=np.array(node_map).T,#nodes of element vertices
                cell_type = [cell_type],#according to vtk format
                order_nodes = order_nodes)
    
    mesh.addAttribute(zone,'zones')
    
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
    npere = int(header.split()[1])
    
    
    node_map = np.zeros((numel,npere),dtype=int)
    #np.array([[0]*numel]*npere,dtype=int) # connection matrix mapping elements onto nodes 
    elm_no = [0]*numel # element number / index 
    zone = [0]*numel # mesh zone 
    
    for i in range(numel):
        line = fh.readline().split()#read in line data
        elm_no[i] = int(line[0])
        zone[i] = int(line[-1])
        
        for j in range(npere):
            node_map[i,j]=int(line[j+1])-1
    ct = 10        
    if npere==3:
        ct = 5

    #create mesh instance 
    mesh = Mesh(node_x = node_x,#x coordinates of nodes 
                node_y = node_y,#y coordinates of nodes
                node_z = node_z,#z coordinates of nodes 
                node_data=node_map,#nodes of element vertices
                cell_type = [ct],#according to vtk format
                original_file_path = file_path,
                order_nodes = order_nodes)
    
    mesh.addAttribute(zone,'zone')
    
    return mesh
        
        
#%% build a quad mesh        
def quadMesh(elec_x, elec_z, elec_type = None, elemx=4, xgf=1.5, zf=1.1, zgf=1.25, fmd=None, pad=2, 
              surface_x=None,surface_z=None,refine_x = None, refine_z=None):
    """Creates a quaderlateral mesh given the electrode x and y positions.
            
    Parameters
    ----------
    elec_x : list, np array
        Electrode x coordinates 
    elec_z : list, np array
        Electrode y coordinates
    elec_type: list, optional
        strings, where 'electrode' is a surface electrode; 'buried' is a buried electrode
    elemx : int, optional
        Number of nodes between electrodes. 
    xgf : float, optional 
        X factor multiplier for fine zone. 
    zf : float, optional
         Z factor multiplier in the fine zone (must be >1).
    zgf : float, optional
         Z factor multiplier in the coarse zone (must be >1).
    fmd : float (m), optional 
         Fine mesh region depth specifies as positive number (if None, half survey width is used).
    pad : int, optional
         X padding outside the fine area (tipicaly twice the number of elements between electrodes).
    surface_x: array like, optional
        Default is None. x coordinates of extra surface topography points for the generation of topography in the quad mesh
    surface_z: array like, optional
        Default is None. z coordinates of extra surface topography points for the generation of topography in the quad mesh. Note
        an error will be returned if len(surface_x) != len(surface_z)
    refine_x: array like, optional
        Default is None. Inserts points in the resulting quad mesh for more control over 
        mesh refinement at depth. X coordinates. An error will be returned if 
        len(refine_x) != len(refine_z).
    refine_z: array like, optional
        Default is None. Inserts points in the resulting quad mesh for more control over 
        mesh refinement at depth. Z coordinates. An error will be returned if 
        len(refine_x) != len(refine_z).
        
            
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
        
    if refine_x is None or refine_z is None: # flag if to use refinement points  
        refine_flag = False 
    else:
        refine_flag = True 
    
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
        norm_bhy = np.interp(bh[:,0], X, Y) - bh[:,1]
        meshz = np.unique(np.append(meshz,norm_bhy))
        
    if refine_flag:
        norm_rz = np.interp(refine_x, X, Y) - refine_z 
        meshz = np.unique(np.append(meshz,norm_rz))
    
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
                np.array(node_mappins).T,
                cell_type=[9],
                original_file_path='N/A',
                order_nodes=False) 
 
    
    #find the node which the electrodes are actually on in terms of the mesh. 
    node_in_mesh = [0]*len(elec_x)
    for i in range(len(elec_x)):
        sq_dist = (node_x - elec_x[i])**2 + (node_z - elec_z[i])**2 # find the minimum square distance
        
        node_in_mesh[i] = np.argmin(sq_dist) # min distance should be zero, ie. the node index.
    mesh.setElecNode(node_in_mesh) # add nodes to the mesh class

    # OVERWRITE elec_node using node_in_mesh that doesn't assume surface
    # electrode but snap them to node based on distance
    elec_node = np.c_[np.arange(len(node_in_mesh))+1, node_in_mesh]

    return mesh, meshx, meshz, topo, elec_node

def quad_mesh(*args):
    warnings.warn('quad_mesh is depreciated, use quadMesh instead')
    quadMesh(*args)


#%% handling gmsh
def runGmsh(ewd, file_name, show_output=True, dump=print, threed=False, handle=None):
    """

    Parameters
    ----------
    ewd : str
        Directory where gmsh copy is stored.
    file_name : str
        Name of the .geo file without extension.
    show_output : TYPE, optional
        If True, output of gmsh is displayed to dump. The default is True.
    threed : bool, optional
        If True, 3D mesh is done, else 2D. The default is False.
    handle : variable, optional
        Will be assigned the output of 'Popen' in case the process needs to be
        killed in the UI for instance.

    Returns
    -------
    None.

    """
    opt = '-2' 
    if threed: # if 3d use 3d option 
        opt = '-3'

    if platform.system() == "Windows":#command line input will vary slighty by system 
        cmd_line = [os.path.join(ewd,'gmsh.exe'), file_name+'.geo', opt, 'nt %i'%ncores]
        
    elif platform.system() == 'Darwin': # its a macOS 
        winetxt = 'wine'
        if getMacOSVersion():
            winetxt = 'wine64'
        winePath = []
        wine_path = Popen(['which', winetxt], stdout=PIPE, shell=False, universal_newlines=True)#.communicate()[0]
        for stdout_line in iter(wine_path.stdout.readline, ''):
            winePath.append(stdout_line)
        if winePath != []:
            cmd_line = ['%s' % (winePath[0].strip('\n')), ewd+'/gmsh.exe', file_name+'.geo', opt,'-nt','%i'%ncores]
        else:
            try:
                is_wine = Popen(['/usr/local/bin/%s' % winetxt,'--version'], stdout=PIPE, shell = False, universal_newlines=True)
                wPath = '/usr/local/bin/'
            except:
                is_wine = Popen(['/opt/homebrew/bin/%s' % winetxt,'--version'], stdout=PIPE, shell = False, universal_newlines=True) # quick fix for M1 Macs
                wPath = '/opt/homebrew/bin/'
            cmd_line = [wPath + winetxt, ewd+'/gmsh.exe', file_name+'.geo', opt,'-nt','%i'%ncores]
    
    elif platform.system() == 'Linux':
        if os.path.isfile(os.path.join(ewd,'gmsh_linux')):
            cmd_line = [ewd + '/gmsh_linux', file_name + '.geo', opt,'-nt','%i'%ncores] # using linux version if avialable (can be more performant)
        else: # fall back to wine
            cmd_line = ['wine',ewd+'/gmsh.exe', file_name+'.geo', opt,'-nt','%i'%ncores]
    
    else:
        raise Exception('Unsupported operating system') # if this even possible? BSD maybe. 

    if show_output: 
        p = Popen(cmd_line, stdout=PIPE, stderr=PIPE, shell=False)#run gmsh with ouput displayed in console
        if handle is not None:
            handle(p)
        while p.poll() is None:
            line = p.stdout.readline().rstrip()
            if line.decode('utf-8') != '':
                dump(line.decode('utf-8'))
    else:
        p = Popen(cmd_line, stdout=PIPE, stderr=PIPE, shell=False)
        if handle is not None:
            handle(p)
        p.communicate() # wait to finish
        
#%% handle repeated nodes 
def check4repeatNodes(X,Y,Z,flag=None):
    """Raise error if repeated nodes present 

    Parameters
    ----------
    X : array like 
        X coordinates of nodes (normally electrodes).
    Y : array like 
        Y coordinates of nodes (normally electrodes).
    Z : array like 
        Z coordinates of nodes (normally electrodes).
    flag : list of string, optional
        String flag assigned to each node. The default is None.

    Raises
    ------
    ValueError
        If repeated nodes detected 

    Returns
    -------
    None.

    """
    if len(X) != len(Y) or len(X) != len(Z):
        raise ValueError('XYZ Arrays not of equal length')
        
    if flag is None:
        flag = ['node']*len(X)
        
    pts = np.c_[X,Y,Z]
    tree = cKDTree(pts)
    dist,neigh = tree.query(pts,2)
    if any(dist[:,1]<1e-15): # raise an issue when repeated nodes are present 
        pblm_idx = np.arange(pts.shape[0])[dist[:,1] < 1e-15]
        
        ## go through pairs of repeated nodes and keep one of them, so pair up the repeated nodes 
        pairs = [] # stores pairs of repeated points 
        probed = [] # chaches any points already probed (hence can be ignored)
        del_idx = [] # indexes of points to be deleted (?)
        for i in pblm_idx:
            pair = (i,neigh[i,1])
            c = 0  # triggers if pairing not already in list 
            if pair[0] not in probed:
                probed.append(pair[0])
                c += 1 
            if pair[1] not in probed:
                probed.append(pair[1])
                c += 1
            if c!=0: 
                pairs.append(pair)

        ## print out to user the problematic nodes 
        error = 'The following surface nodes are repeated!\n'
        # print(error.strip())
        for pair in pairs:
            for p in pair: 
                x_ = X[p]
                y_ = Y[p]
                z_ = Z[p]
                f_ = flag[p]
                line = "X = {:16.8f}, Y = {:16.8f}, Z = {:16.8f}, flag = {:}\n".format(x_,y_,z_,f_)
                # print(line.strip())
                error += line 
        
        #warnings.warn(error,Warning)
        raise ValueError(error)


#%% build a triangle mesh - using the gmsh wrapper
def triMesh(elec_x, elec_z, elec_type=None, geom_input=None, keep_files=True, 
             show_output=True, path='exe', dump=print, whole_space=False, 
             handle=None, **kwargs):
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
        - 'remote' = remote electrode (not actually placed in the mesh)
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
    handle : variable, optional
        Will be assigned the output of 'Popen' in case the process needs to be
        killed in the UI for instance.
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
    
    # run gmsh
    runGmsh(ewd, file_name, show_output=show_output, dump=dump, threed=False, handle=handle)
        
    #convert into mesh.dat
    mesh_info = gw.mshParse(file_name+'.msh', debug=show_output) # read in mesh file
    
    # merge fine with coarse regions (coarse = 1, fine = -1)
    regions = np.array(mesh_info['parameters'])
    ureg = np.sort(np.unique(regions))
    ie = regions == ureg[-1] # fine region
    regions[ie] = 1 # is same number as the coarse
    
    mesh = Mesh(mesh_info['node_x'], # convert output of parser into an object
                mesh_info['node_y'],
                mesh_info['node_z'],
                np.array(mesh_info['node_data']).T,
                mesh_info['cell_type'],
                mesh_info['original_file_path'])
    
    mesh.addAttribute(regions, 'region')
    
    if keep_files is False: 
        os.remove(file_name+".geo")
        os.remove(file_name+".msh")

    mesh.setElecNode(node_pos-1)#in python indexing starts at 0, in gmsh it starts at 1 
    
    if elec_type is not None:
        iremote = np.array([a == 'remote' for a in elec_type])
        mesh.iremote = iremote

    return mesh

def tri_mesh(*args):
    warnings.warn('tri_mesh is depreciated, use triMesh instead')
    triMesh(**args)



#%% 3D tetrahedral mesh 
def tetraMesh(elec_x,elec_y,elec_z=None, elec_type = None, keep_files=True, interp_method = 'triangulate',
               surface_refinement=None, mesh_refinement=None, show_output=True, 
               path='exe', dump=print, whole_space=False, padding=20,
               search_radius = 10, handle=None, **kwargs):
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
    show_output : boolean, optional
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
    handle : variable, optional
        Will be assigned the output of 'Popen' in case the process needs to be
        killed in the UI for instance.
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
        
    # check for repeated electrodes? 
    check4repeatNodes(elec_x,elec_y,elec_z,elec_type)
            
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
            dump('found buried electrodes')
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
         
    # check directories 
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
        dump("Whole space problem")
        raise Exception("Sorry whole space 3D problems are not implemented yet")
        
    else:
        node_pos = gw.box_3d([elec_x,elec_y,elec_z], file_path=file_name, **kwargs)
            
    # handling gmsh
    runGmsh(ewd, file_name, show_output=show_output, dump=dump, threed=True, handle=handle)
        
    #convert into mesh.dat
    mesh_info = gw.mshParse(file_name+'.msh', debug=show_output) # read in 3D mesh file
    
    # merge fine with coarse regions
    regions = np.array(mesh_info['parameters'])
    for reg in np.unique(regions)[1:]:
        ie = regions == reg
        regions[ie] = reg - 1
    
    mesh = Mesh(mesh_info['node_x'], # convert output of parser into an object
                mesh_info['node_y'],
                mesh_info['node_z'],
                np.array(mesh_info['node_data']).T,
                mesh_info['cell_type'],
                mesh_info['original_file_path'])

    mesh.addAttribute(regions, 'region')
    mesh.addAttribute(np.array(mesh_info['parameters']), 'gmsh_phys_entity')
    
    #mesh.write_dat(file_path='mesh.dat') # write mesh.dat - disabled as handled higher up in the R2 class 
    node_x = np.array(mesh.node[:,0])
    node_y = np.array(mesh.node[:,1])
    
    if keep_files is False: 
        os.remove(file_name+".geo");os.remove(file_name+".msh")
        
    dump('interpolating topography onto mesh using %s interpolation...'%interp_method)
    
    x_interp = np.append(surf_elec_x,surf_x)#parameters to be interpolated with
    y_interp = np.append(surf_elec_y,surf_y)
    z_interp = np.append(surf_elec_z,surf_z)

    #using home grown functions to interpolate / extrapolate topography on mesh
    if interp_method == 'idw': 
        nodez = interp.idw(node_x, node_y, x_interp, y_interp, z_interp,radius=search_radius)# use inverse distance weighting
    elif interp_method == 'bilinear':# interpolate on a irregular grid, extrapolates the unknown coordinates
        nodez = interp.interp2d(node_x, node_y, x_interp, y_interp, z_interp)
    elif interp_method == 'nearest':
        nodez = interp.nearest(node_x, node_y, x_interp, y_interp, z_interp,num_threads=ncores)
    elif interp_method == 'spline':
        nodez = interp.interp2d(node_x, node_y, x_interp, y_interp, z_interp,method='spline')
    elif interp_method == 'triangulate':
        nodez = interp.triangulate(node_x, node_y, x_interp, y_interp, z_interp)
    elif interp_method == None:
        nodez = np.zeros_like(node_x,dtype=float)
    
    mesh.node[:,2] = mesh.node[:,2] + nodez
    #need to recompute cell centres as well as they will have changed. 
    mesh.cellCentres()

    #check if remeote electrodes present, and insert them into the node position array
    if len(rem_elec_idx)>0: 
        rem_node_bool = min(node_x) == node_x #& (min(node_y) == node_y) & (min(node_z) == node_z)
        rem_node = np.argwhere(rem_node_bool == True)[0][0]
        #node_pos = np.insert(node_pos,np.array(rem_elec_idx),[rem_node]*len(rem_elec_idx),axis=0)
        # insert node position in to electrode node array 
        node_pos_tmp = [] # temporary node position array 
        iremote = [False]*len(elec_type)
        c = 0 
        for i, key in enumerate(elec_type):
            if key == 'remote':
                node_pos_tmp.append(rem_node+1)
                iremote[i] = True 
            else: 
                node_pos_tmp.append(node_pos[c])
                c+=1            
        node_pos = np.array(node_pos_tmp,dtype=int) # redefine node positioning 
        mesh.setElecNode(node_pos-1,np.array(iremote,dtype=bool))
    else:
        #add nodes to mesh
        mesh.setElecNode(node_pos-1)#in python indexing starts at 0, in gmsh it starts at 1 
    
    return mesh

def tetra_mesh(*args):
    warnings.warn('tetra_mesh is depreciated, use tetraMesh instead', DeprecationWarning)
    tetraMesh(*args)



#%% column mesh 
def prismMesh(elec_x, elec_y, elec_z, 
               keep_files=True,
               show_output=True, 
               path='exe', dump=print,
               handle=None, **kwargs):
    """Make a prism mesh.
    
    Parameters
    ----------
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
    handle : variable, optional
        Will be assigned the output of 'Popen' in case the process needs to be
        killed in the UI for instance.
    """
    #error checks 
    if len(elec_x) != len(elec_y):
        raise ValueError("mismatch in electrode x and y vector length, they should be equal")
        
    #make prism mesh command script
    file_name = "prism_mesh"
    gw.prism_mesh([elec_x,elec_y,elec_z], file_path=file_name, **kwargs)
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
    runGmsh(ewd, file_name, show_output=show_output, dump=dump, threed=True, handle=handle)
        
    #convert into mesh.dat
    mesh_info = gw.mshParse(file_name+'.msh', debug=show_output) # read in 3D mesh file
   
    # merge fine with coarse regions
    regions = np.array(mesh_info['parameters'])
    for reg in np.unique(regions)[1:]:
        ie = regions == reg
        regions[ie] = reg - 1
        
    mesh = Mesh(mesh_info['node_x'], # convert output of parser into an object
                mesh_info['node_y'],
                mesh_info['node_z'],
                np.array(mesh_info['node_data']).T,
                mesh_info['cell_type'],
                mesh_info['original_file_path'])
    
    mesh.addAttribute(regions, 'region')
    mesh.moveElecNodes(elec_x,elec_y,elec_z)
    
    if keep_files is False: 
        os.remove(file_name+".geo");os.remove(file_name+".msh")
        
    return mesh 

def prism_mesh(*args):
    warnings.warn('prism_mesh is depreciated, use prismMesh instead', DeprecationWarning)
    prismMesh(*args)
    

#%% cylinder mesh
def cylinderMesh(elec, zlim=None, radius=None, file_path='cylinder_mesh.geo',
               cl=-1, finer=4,
               keep_files=True,
               show_output=True, 
               path='exe', dump=print,
               handle=None):
    """Make a cylinder mesh.
    
    Parameters
    ----------
    elec : list of array_like or nx3 array
        First column/list is the x coordinates of electrode positions and so on ... 
    zlim : list, tuple, optional 
        Bottom and top z coordinate of column, in the form (min(z),max(z))
    radius: float, optional 
        Radius of column. If not provided, will be infered from elec positions.
    file_path : string, optional 
        Name of the generated gmsh file (can include file path also).
    cl : float, optional
        Characteristic length, essentially describes how big the nodes 
        associated elements will be. Usually no bigger than 5. If set as -1 (default)
        a characteristic length 1/4 the minimum electrode spacing is computed. 
    finer : int, optional
        Number of line between two consecutive electrodes to approximate the
        circle shape.
    handle : variable, optional
        Will be assigned the output of 'Popen' in case the process needs to be
        killed in the UI for instance.
    """
    file_path = file_path if file_path[-4:] == '.geo' else file_path + '.geo'
    gw.cylinder_mesh(elec, zlim=zlim, radius=radius, file_path=file_path,
                     cl=cl, finer=finer)    
    
    # check directories 
    if path == "exe":
        ewd = os.path.join(
                os.path.dirname(os.path.realpath(__file__)),
                path)
    else:
        ewd = path # points to the location of the .exe 
        # else its assumed a custom directory has been given to the gmsh.exe 
    
    # handling gmsh
    runGmsh(ewd, file_path.replace('.geo',''), show_output=show_output, dump=dump, threed=True, handle=handle)        

    #convert into mesh.dat
    mesh_info = gw.mshParse(file_path.replace('.geo', '.msh'), debug=show_output) # read in 3D mesh file

    # merge fine with coarse regions
    regions = np.array(mesh_info['parameters'])
    for reg in np.unique(regions)[1:]:
        ie = regions == reg
        regions[ie] = reg - 1
    regions[regions == 0] = 1 # fix to be sure the first region is 1 not 0
        
    # create mesh object        
    mesh = Mesh(mesh_info['node_x'], # convert output of parser into an object
                mesh_info['node_y'],
                mesh_info['node_z'],
                np.array(mesh_info['node_data']).T,
                mesh_info['cell_type'],
                mesh_info['original_file_path'])
    
    mesh.addAttribute(regions, 'region')
    mesh.moveElecNodes(elec[:,0], elec[:,1], elec[:,2])
    
    if keep_files is False: 
        os.remove(file_path)
        os.remove(file_path.replace('.geo', '.msh'))
        
    return mesh

# test
# radius = 6.5/2 # cm
# angles = np.linspace(0, 2*np.pi, 13)[:-1] # radian
# celec = np.c_[radius*np.cos(angles), radius*np.sin(angles)]
# elec = np.c_[np.tile(celec.T, 8).T, np.repeat(6.5+np.arange(0, 8*5.55, 5.55)[::-1], 12)]
# mesh = cylinderMesh(elec, file_path='/home/jkl/Downloads/mesh_cylinder',
#                     zlim=[0, 47.5], cl=0.5)

#%% tank mesh
def tankMesh(elec=None, 
             origin=None, 
             dimension=[10,10,10],
             file_path='tank_mesh.geo',
             cl=-1,
             keep_files=True,
             show_output=True, 
             path='exe',
             dump=print,
             handle=None):
    """Make a tank mesh (3D closed box).
    
    Parameters
    ----------
    elec : list of array_like or nx3 array
        First column/list is the x coordinates of electrode positions and so on ... 
    origin : list of float, optional
        Origin of the corner where the mesh will be drawned from. If not
        provided and elec provided, the smaller elec position will be chosen.
    dimension : list of float, optional
        Dimension of the mesh in X,Y,Z from the corner origin.
        The default is [10,10,10].
    file_path : string, optional 
        Name of the generated gmsh file (can include file path also).
    cl : float, optional
        Characteristic length, essentially describes how big the nodes 
        associated elements will be. Usually no bigger than 5. If set as -1 (default)
        a characteristic length 1/4 the minimum electrode spacing is computed. 
    finer : int, optional
        Number of line between two consecutive electrodes to approximate the
        circle shape.
    handle : variable, optional
        Will be assigned the output of 'Popen' in case the process needs to be
        killed in the UI for instance.
    """
    file_path = file_path if file_path[-4:] == '.geo' else file_path + '.geo'
    gw.tank_mesh(elec=elec, origin=origin, dimension=dimension,
                 file_path=file_path, cl=cl)    
    
    # check directories 
    if path == "exe":
        ewd = os.path.join(os.path.dirname(os.path.realpath(__file__)), path)
    else:
        ewd = path # points to the location of the .exe 
        # else its assumed a custom directory has been given to the gmsh.exe 
    
    # handling gmsh
    runGmsh(ewd, file_path.replace('.geo',''), show_output=show_output, dump=dump, threed=True, handle=handle)        

    #convert into mesh.dat
    mesh_info = gw.mshParse(file_path.replace('.geo', '.msh'), debug=show_output) # read in 3D mesh file

    # merge fine with coarse regions
    regions = np.array(mesh_info['parameters'])
    for reg in np.unique(regions)[1:]:
        ie = regions == reg
        regions[ie] = reg - 1
        
    # create mesh object        
    mesh = Mesh(mesh_info['node_x'], # convert output of parser into an object
                mesh_info['node_y'],
                mesh_info['node_z'],
                np.array(mesh_info['node_data']).T,
                mesh_info['cell_type'],
                mesh_info['original_file_path'])
    
    mesh.addAttribute(regions, 'region')
    mesh.moveElecNodes(elec[:,0], elec[:,1], elec[:,2])
    
    if keep_files is False: 
        os.remove(file_path)
        os.remove(file_path.replace('.geo', '.msh'))
        
    return mesh

# test
# TODO debug
# elec = np.array([[0,2,2],[0,2,6],[0,3,2],[0,3,6],
#                   [10,2,2],[10,2,6],[10,3,2],[10,3,6],
#                   [3,0,2],[5,0,2],[7,0,2],[3,0,6],[5,0,6],[7,0,6],
#                   [3,5,2],[5,5,2],[7,5,2],[3,5,6],[5,5,6],[7,5,6]
#                   ])
# mesh = tankMesh(elec, origin=[0,0,0], dimension=[10, 5, 7], file_path='/home/jkl/Downloads/mesh_tank',)
# mesh.show()

#%% ciruclar mesh (TODO)
def circularMesh():
    # TODO
    pass

#%% import a custom mesh, you must know the node positions 
def readMesh(file_path, node_pos=None, order_nodes=True):
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
        mesh_dict = gw.mshParse(file_path, debug=False)

        mesh = Mesh(node_x = mesh_dict['node_x'],
                    node_y = mesh_dict['node_y'],
                    node_z = mesh_dict['node_z'],
                    node_data = np.array(mesh_dict['node_data']).T,
                    cell_type = mesh_dict['cell_type'],
                    order_nodes = order_nodes)
        mesh.addAttribute(mesh_dict['parameters'], 'region')
        
    elif ext == '.dat':
        mesh = dat_import(file_path, order_nodes=order_nodes)   
    elif ext == '.node':
        mesh = tetgen_import(file_path, order_nodes=order_nodes)
    else:
        avail_ext = ['.vtk','.msh','.dat','.node']
        raise ImportError("Unrecognised file extension, available extensions are "+str(avail_ext))
    
    if node_pos is not None:
        mesh.setElecNode(np.array(node_pos, dtype=int)) # add electrode nodes to mesh provided by the user
    
    return mesh

def custom_mesh_import(file_path, node_pos=None, order_nodes=True):
    warnings.warn('custom_mesh_import is depreciated, use readMesh instead', DeprecationWarning)
    mesh = readMesh(file_path, node_pos, order_nodes)
    return mesh
    
#%% merge meshes from multiple mesh objects 
def mergeMeshes(meshList):
    """Merge multiple meshes into one mesh object, useful for psuedo 3d surveys. 

    Parameters
    ----------
    meshList : list
        List of seperate mesh classes .

    Returns
    -------
    mmesh: class
        Merged mesh. 

    """
    mesh = meshList[0]
    conmat = mesh.connection.copy()
    node = mesh.node.copy()
    centriod = mesh.elmCentre.copy()
    df = mesh.df.copy()
    ct = np.array(mesh.cell_type)
    iremote = mesh.iremote
    zone = mesh.zone 
    numnp = np.max(conmat.flatten())
    
    #electrode handling
    try:
        if mesh.eNodes is not None:
            eNodes = mesh.eNodes.copy()
        else:
            eNodes = None
    except:
        eNodes = None
    
    iremote_present = True
    if iremote is None:
        iremote_present = False 
    zone_present = True
    if zone is None:
        zone_present = False 
    node_present = True
    if eNodes is None:
        node_present = False
    
    for i in range(1,len(meshList)):
        mesh = meshList[i]
        nconmat = mesh.connection + numnp + 1
        conmat = np.vstack([conmat,nconmat])
        node = np.vstack([node,mesh.node])
        centriod = np.vstack([centriod, mesh.elmCentre])
        df = pd.concat([df,mesh.df],ignore_index=True)
        ct = np.concatenate([ct,np.array(mesh.cell_type)],axis=0)
        if iremote_present:
            iremote = np.concatenate([iremote,mesh.iremote])
        if zone_present:
            zone = np.concatenate([zone,mesh.zone]) 
        if node_present:
            eNodes = np.concatenate([eNodes,mesh.eNodes+numnp+1])
        numnp = np.max(conmat.flatten())
        
    mmesh = Mesh(node[:,0],
                 node[:,1],
                 node[:,2],
                 conmat,
                 ct,
                 'merged from several meshes',
                 order_nodes=False,# should be already ordered
                 compute_centre=False,# should be already computed 
                 check2D=False)# should be already computed 
    
    mmesh.elmCentre = centriod
    mmesh.df = df.copy()
    
    mmesh.iremote = iremote
    mmesh.zone = zone
        
    return mmesh # merged mesh  