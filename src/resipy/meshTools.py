# MESH TOOLS 
"""
Created on Wed May 30 10:19:09 2018, python 3.6.5
@author: jamyd91
Module handles mesh generation, display, discretisation and post processing. 
The convention for x y z coordinates is that the z coordinate is the elevation.

Dependencies: 
    numpy (standard scientific lib)
    matplotlib (standard scientific lib)
    gmshWrap(ResIPy resipy module)
    python3 standard libaries
"""
#import standard python packages
import os, platform, warnings, psutil, struct, base64, time, ntpath, shutil
from subprocess import PIPE, Popen, run, TimeoutExpired
import tempfile
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import ListedColormap
import matplotlib.tri as tri
import matplotlib.patches as mpatches
import matplotlib.path as mpath
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import cKDTree, Voronoi, ConvexHull, Delaunay
from scipy.interpolate import interp1d
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

def findDist(elec_x, elec_y, elec_z): # find maximum and minimum electrode spacings 
    dist = np.zeros((len(elec_x),len(elec_x)))   
    x1 = np.array(elec_x)
    y1 = np.array(elec_y)
    z1 = np.array(elec_z)
    for i in range(len(elec_x)):
        x2 = elec_x[i]
        y2 = elec_y[i]
        z2 = elec_z[i]
        dist[:,i] = np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
    return dist.flatten() # array of all electrode distances 

def findminmax(a, pad=20):
    a = np.array(a)
    delta = abs(np.max(a) - np.min(a))
    mina = np.min(a) - (delta*(pad/100))
    maxa = np.max(a) + (delta*(pad/100))
    
    if mina == maxa:
        mina -= 1
        maxa += 1
    return [mina, maxa]

#%% check mac version for wine
def getMacOSVersion(): # this is obsolete now
    OpSys=platform.system()    
    if OpSys=='Darwin':
        versionList = platform.mac_ver()[0].split('.')
        macVersion = float(versionList[0] + '.' + versionList[1]) # not getting patch version so xx.xx only
        if macVersion >= 10.15:
            return True
        else:
            return False

def whichWineMac():
    """
    Checks if 'wine' or 'wine64' should be used on macOS.

    Returns:
        str: The appropriate wine executable ('wine' or 'wine64'), or None if neither is found.
    """    

    wine_paths = ['/usr/local/bin/wine', '/opt/homebrew/bin/wine', '/usr/bin/wine'] # common wine paths
    wine64_paths = ['/usr/local/bin/wine64', '/opt/homebrew/bin/wine64', '/usr/bin/wine64'] # common wine64 paths
    
    global wPath
    for path in wine_paths:
        if os.path.exists(path) and os.access(path, os.X_OK):
            try:
                # Basic check if wine is working
                run([path, '--version'], capture_output=True, timeout=1)
                wPath = path
                return 'wine'
            except (TimeoutExpired, FileNotFoundError, OSError):
                pass # wine either timed out, was not found, or had an OS error.

    for path in wine64_paths:
        if os.path.exists(path) and os.access(path, os.X_OK):
            try:
                # Basic check if wine64 is working
                run([path, '--version'], capture_output=True, timeout=1)
                wPath = path
                return 'wine64'
            except (TimeoutExpired, FileNotFoundError, OSError):
                pass # wine64 either timed out, was not found, or had an OS error.

    return None
        
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
    -------
    Mesh : class
    
    Attritbutes 
    -----------
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
    
    Note
    ----
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
        
        self.dint = np.int64 # connection matrix needs to be in long format 
        if platform.system() == 'Windows':
            self.dint = np.int32 # avoid windows quirk where type long is actually a 32 bit integer
        elif platform.machine() == 'armv7l':
            self.dint = np.int64
        self.connection = np.asarray(node_data,dtype=self.dint) #connection matrix
        
        self.cell_type = cell_type # cellType
        self.originalFilePath = original_file_path # originalFilePath 
        self.eNodes  = None # node number that electrodes occupy 
        self.elec = None
        self.elec_type = None 
        self.fmd = None 
        self.iremote = None # specify which electrode is remote
        self.idirichlet = None # specify idirichlet node 
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
        
        # number of voxels in the xyz directions, only relevant to the voxel mesh type 
        self.nvoxel = {'x':None, 'y':None, 'z':None}
    
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
        ----------
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
        elif int(self.cell_type[0]) == 11 or int(self.cell_type[0]) == 12: # elements are voxels
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
        elif int(self.cell_type[0]) == 11 or int(self.cell_type[0]) == 12:# elements are voxels
            return 8
        elif int(self.cell_type[0]) == 10:# elements are tetrahedra 
            return 4
        elif int(self.cell_type[0]) == 13: # elements are 3d wedges 
            return 5
        #add element types as neccessary 
        else:
            warnings.warn("Unrecognised cell type")
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
        ----------
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
        self.df[key] = values  #allows us to add an attributes to each element.
        
    def addPtAttribute(self,values,key):
        """Associate attributes with the mesh nodes  
        
        Parameters
        ----------
        values: array like
            Must have a length which matches the number of nodes. Discrete 
            values which map to the nodes. 
        key: str
            Name of the attribute, this will be used to reference to the values
            in other mesh functions. 
        """
        if len(values)!=self.numnp:
            raise ValueError("The length of the new attributes array (%i) does not match the number of nodes in the mesh (%i)"%(len(values),self.numnp))
        if key in self.ptdf.keys():
            self.ptdf[key] = values # reallocate column if already set 
        else:
            self.ptdf[key]=values #allows us to add an attributes to each point 
            
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
        ----------
        return_count:bool, optional
            Function returns the number of reordered elements, default is False. 
        """
        con_mat = self.connection
        self.dint = self.connection.dtype
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

        self.connection = np.asarray(con_mat,dtype=self.dint)
        
        if return_count:
            return count
        
    def orderElem(self, param=None):
        """Order elements based on the parameter number. Ideally parameters 
        should be concurrent, and fixed elements should go at the base of the 
        connection matrix.
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
        self.computeNeigh()
        
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
                self.neigh_matrix, self.tri_combo = mc.neigh3d(np.asarray(self.connection, dtype=self.dint),1,ncores)
            elif self.type2VertsNo() == 6:  # prism mesh 
                self.neigh_matrix, self.tri_combo = mc.neighPrism(self.connection,1,ncores)


    def computeNconnec(self):
        cell_type = self.cell_type[0]
        self.NsizeA, self.fconm = mc.conductanceCall(self.connection, self.numnp,
                                                     cell_type,ncores)
        
    def externalNodes(self):
        """
        Find the external nodes of the mesh (triangle and tetrahedral meshes only). 

        Returns
        -------
        external_nodes : ndarray 
            Indices of external nodes. 
        surface_flag : ndarray 
            Flag if external node is on the surface of the mesh, 1 for surface
            and 0 for not. 

        """
        throw_error = True 
        if self.ndims==2 and self.type2VertsNo()==3:
            throw_error = False 
        elif self.ndims==3 and self.type2VertsNo()==4:
            throw_error = False 
        if throw_error:
            raise Exception('External boundary node call only avialable for triangle or tetrahedral meshes')
        if self.neigh_matrix is None: 
            self.computeNeigh()
        external_nodes, surface_flag = mc.externalN(self.connection, self.node,
                                                    self.neigh_matrix) 
        return external_nodes, surface_flag
        
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
            print('Sorry, refine method not implemented for this mesh type')
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
        ----------
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
            typ = 2 # flag the type of algorithm to use 
            if self.type2VertsNo() == 4: 
                qmesh = self.quad2tri() # get the triangle mesh instead 
                return qmesh.extractSurface() # run the same function 
        elif self.ndims == 3: 
            typ = 3
            if self.type2VertsNo() == 6 or self.type2VertsNo() == 8:#not a tetrahedra 3d mesh 
                typ = 4 
        
        if typ != 4: # below parameters not used for type 4 flag 
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
            
        elif typ == 3:  ### Extract 3D faces of top of the mesh ###
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
            
        else: 
            # in this case assume that the nodes are all aligned in the horizontal axis 
            # firstly find the unique x y positions in the mesh 
            # from https://stackoverflow.com/questions/7989722/finding-unique-points-in-numpy-array
            unixy = np.vstack([np.array(u) for u in set([tuple(p) for p in self.node[:,0:2]])])
            unode = np.zeros((unixy.shape[0],3), dtype=float)

            # secondly triangulate the nodes at the uppermost position of each 
            # vertical column 
            for i in range(unode.shape[0]):
                x = unixy[i,0]
                y = unixy[i,1]
                boolidx = (x==self.node[:,0]) & (y==self.node[:,1])
                zvals = self.node[boolidx,2] # values of z in the xy column 
                z = max(zvals) # take the maximum value as the upper most point 
                unode[i,0] = x; unode[i,1] = y; unode[i,2] = z; 
                
            # thirdly triangulate between the unique nodes in the xy plane 
            tri = Delaunay(unode[:,0:2]) # triangulate 
            
            nmesh = Mesh(unode[:,0], # make new mesh 
                         unode[:,1], 
                         unode[:,2], 
                         node_data = tri.simplices, 
                         cell_type = [5], 
                         order_nodes=False,
                         compute_centre=True,
                         check2D=False)
            
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
        if len(sortf) == 0:
            warnings.warn('call made to remove excess nodes but there are no excess nodes to remove')
            return 
        
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
        """Crop the mesh given a polyline in 2D. If 3D then the mesh will be 
        cropped in the XY plane. 
        
        Parameters
        ----------
        polyline : array of float
            Array of size Nx2 with the XZ (or XY) coordinates forming the 
            polyline. Notethat the first and last coordinates should be the
            same to close the polyline.
        """
        # get points inside the polygon
        path = mpath.Path(polyline)
        
        if self.ndims == 2: 
            centroids = np.c_[self.elmCentre[:,0], self.elmCentre[:,2]]
        else:
            centroids = np.c_[self.elmCentre[:,0], self.elmCentre[:,1]]
            
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
        ----------
        xlim : tuple, optional
            Axis x limits as `(xmin, xmax)`.
        ylim : tuple, optional
            Axis y limits as `(ymin, ymax)`. 
        zlim : tuple, optional
            Axis z limits as `(ymin, ymax)`. 
            
        Returns
        -------
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
            
        # catch 2d truncation in 3D meshes 
        if ylim[0] == ylim[1]:
            ylim[0] -= 3.0
            ylim[1] += 3.0
        if xlim[0] == xlim[1]:
            xlim[0] -= 3.0
            xlim[1] += 3.0
        if zlim[0] == zlim[1]:
            zlim[0] -= 3.0
            zlim[1] += 3.0
             
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
        ----------
        attr: string
            Name of attribute to threshold by 
        vmin: float
            minimum value of attribute
        vmax: float
            maximum value of attribute
        
        Returns
        -------
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
        ----------
        in_elem: array like
            array of bool, True means to keep the element 
        attr: string
            Name of attribute to threshold by 
        
        Returns
        -------
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
        
        nmesh.df = new_df.reset_index()
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
        ----------
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
        
        if '(log10)' in attr: # the real values (later in the code) cause confusion
            if 'Resistivity' in attr:
                color_bar_title = attr.replace('(log10)', '(ohm.m)')
            else:
                color_bar_title = attr.replace('(log10)', '')
        else:
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
        self.fig.set_tight_layout(True) # set figure to tight layout 
        
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
                # levels = MaxNLocator().tick_values(vmin, vmax)
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
            self.cbar = plt.colorbar(self.cax, ax=ax, orientation=cbar_horizontal, fraction=0.046, pad=0.04)
            if attr == 'region': # default to material
                val = np.sort(np.unique(X)).astype(int)
                if len(val) > 1:
                    interval = (val[-1]-val[0])/len(val)
                    self.cbar.set_ticks(np.arange(val[0]+interval/2, val[-1], interval))
                else:
                    self.cbar.set_ticks([1])
                self.cbar.set_ticklabels(val)
            elif 'log' in attr:
                # levels = MaxNLocator().tick_values(vmin, vmax)
                levels = np.linspace(vmin, vmax, 13)
                if vmax < 0.01:
                    ticks = ['{:.2e}'.format(10**lvl) for lvl in levels]
                else:
                    ticks = ['{:.2f}'.format(10**lvl) for lvl in levels]
                self.cbar.set_ticks(levels)
                self.cbar.set_ticklabels(ticks)
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
        elif 'log' in attr: 
            # levels = MaxNLocator().tick_values(vmin, vmax)
            levels = np.linspace(vmin, vmax,13)
            if vmax < 0.01:
                ticks = ['{:.2e}'.format(10**lvl) for lvl in levels]
            else:
                ticks = ['{:.2f}'.format(10**lvl) for lvl in levels]
            self.cbar.set_ticks(levels)
            self.cbar.set_ticklabels(ticks)
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
    
    
    def show3D(self,color_map = 'Spectral',#displays the mesh using matplotlib or pyvista 
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
                elec_size=10.,
                use_pyvista=True,
                background_color=(0.8,0.8,0.8),
                pvslices=([],[],[]),
                pvspline=None,
                pvthreshold=None,
                pvgrid=True,
                pvcontour=[],
                pseudo3DContour=False,
                pvdelaunay3d=False,
                pvshow=True,
                volume = None, 
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
        elec_size : float, optional
            Size of electrode points, default is 10, a number of 4-5 gives better visibility
        use_pyvista : bool, optional
            Use visual toolkit backend for displaying 3D mesh, note that pyvista
            must be installed for this to work. 
        background_color : tuple, optional 
            Background color assigned to pyvista plotter object when created. Not yet
            supported for matplotlib axis handles. 
        pvslices : tuple of list of float, optional
            Determine the X, Y, Z slices. e.g.: ([3], [], [-3, -4]) will add
            a slice normal to X in 3 and two slices normal to Z in -3 and -4.
        pvspline : 'elec' or numpy array, optional
            If 'elec' mesh will be sliced along the electrodes path otherwise 
            an array of X, Y, Z of points on a path to slice the mesh along that path is needed.
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
        volume : float, optional
            If not 'None' then volume float number will be printed onto the pyvista plot.
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
            
            # apply threshold (do this before sending to pyvista now)
            if pvthreshold is not None:
                if len(pvthreshold) < 2:
                    raise Exception('pvthreshold argument should be array like with at least two entries')
                if pvthreshold[0] is None:
                    pvthreshold[0] = np.nanmin(X)
                if pvthreshold[1] is None:
                    pvthreshold[1] = np.nanmax(X)
                nmesh = nmesh.threshold(attr=color_bar_title, vmin=pvthreshold[0], vmax=pvthreshold[1])
                # self.pvmesh = self.pvmesh.threshold(value=pvthreshold, scalars=attr)
            
            #save to temperory directory 
            folder = tempfile.TemporaryDirectory()
            fname = os.path.join(folder.name, '__to_pv_mesh.vtu')
            nmesh.vtu(fname)
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
            elif pvspline is not None:
                if pvspline == 'elec':
                    pvspline = self.elec.copy()[~iremote,:]
                spline_sliced_path = pv.Spline(pvspline, int(len(pvspline)*10)) # pyvista spline for slicing (the number is how smooth we want it *10 makes 10 times more points on the slice path - smoother)
                mesh_spline_slice = self.pvmesh.slice_along_line(spline_sliced_path)
                if mesh_spline_slice.number_of_points > 0 or mesh_spline_slice.number_of_cells > 0:
                    ax.add_mesh(mesh_spline_slice,
                                scalars = attr,
                                cmap=color_map,
                                clim=[vmin, vmax],
                                show_scalar_bar=color_bar,
                                show_edges=edges,
                                opacity=alpha,
                                scalar_bar_args={'color':tcolor},
                                lighting=True)
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
                                lighting=True,
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
                    ax.add_mesh(pvelec, color=elec_color, point_size=elec_size,
                                render_points_as_spheres=True)
                except AttributeError as e:
                    print("Could not plot 3d electrodes, error = "+str(e))
            
            # Volume for 3D
            if volume is not None:
                ax.add_text('Volume = {:.2f} m³'.format(volume), position='upper_left', color=tcolor)    
            
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
        elif self.type2VertsNo() == 8:
            raise Exception('Cannot show a voxel mesh using matplotlib calls, in this version of ResIPy')
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
        fcon, idx = mc.faces3d(np.asarray(tmesh.connection,dtype=self.dint), 
                               np.asarray(tmesh.neigh_matrix,dtype=self.dint))
        

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
            
        if self.iremote is not None: 
            iremote = self.iremote
        else:
            nelec = self.elec.shape[0]
            iremote = np.array([False]*nelec,dtype=bool)
        # work out extents of mesh surface 
        xlim = [min(self.node[:,0])-1, max(self.node[:,0])+1]
        if self.elec is not None: 
            xlim = findminmax(self.elec[:,0][~iremote])

        ylim = [min(self.node[:,1])-1, max(self.node[:,1])+1]
        if self.elec is not None: 
            ylim = findminmax(self.elec[:,1][~iremote])
            
        zlim = None 
        if self.fmd is not None and self.elec is not None:
            zlim = [max(self.elec[:,2])-self.fmd, max(self.elec[:,2])]
        tmesh = self.truncateMesh(xlim,ylim,zlim)

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
        if self.iremote is not None: 
            iremote = self.iremote
        else:
            nelec = self.elec.shape[0]
            iremote = np.array([False]*nelec,dtype=bool)
        if xlim is None:
            xlim = [min(self.node[:,0])-1, max(self.node[:,0])+1]
            if self.elec is not None: 
                xlim = findminmax(self.elec[:,0][~iremote])
        if ylim is None:
            ylim = [min(self.node[:,1])-1, max(self.node[:,1])+1]
            if self.elec is not None: 
                ylim = findminmax(self.elec[:,1][~iremote])
        if zlim is None: 
            if self.fmd is not None and self.elec is not None:
                zlim = [max(self.elec[:,2])-self.fmd, max(self.elec[:,2])]
            else:
                zlim = [min(self.node[:,2]), max(self.node[:,2])]
                
        # if 'elm_id' not in self.df.keys():
        self.df['elm_id'] = np.arange(1,self.numel+1,1)
        
        # surface mesh
        if self.surfaceMesh is None: # do this to avoid recomputing surface mesh
            self.surfaceMesh = self._select3Dsurface() # surface of mesh 
        self.surfaceMesh.show3D(ax=ax,color_map='Greys', color_bar=False,
                     edge_color=None, alpha=0.8, darkMode=darkMode)
        
        # pyvista mesh
        folder = tempfile.TemporaryDirectory()
        fname = os.path.join(folder.name, '__to_pv_mesh.vtu')
        self.vtu(fname)
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
        idx = np.unique(np.asarray(clipped.cell_data['elm_id'],dtype=self.dint))-1
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
            
            # compute for nodes too
            nZ = np.interp(self.node[:,0],datum_x,datum_z,
                           left=datum_z[min_idx],
                           right=datum_z[max_idx])
            node_depths = nZ - self.node[:,2]
            self.ptdf['depths'] = node_depths
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
            
            # compute for nodes too
            nZ = interp.triangulate(self.node[:,0], self.node[:,1], 
                                    surmesh.node[:,0], 
                                    surmesh.node[:,1], 
                                    surmesh.node[:,2])
            node_depths = nZ - self.node[:,2]
            self.ptdf['depths'] = node_depths
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
        ideally this is implemented on a mesh which is refined near the surface. If 
        no nodes are assigned to mesh, a mesh.e_nodes variable is created.  
        
        Parameters
        ----------
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
        -------
        node_in_mesh: np array
            Array of mesh node numbers corresponding to the electrode postions/
            coordinates.
            
        Notes
        -----
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
     
        return np.array(node_in_mesh,dtype=self.dint) # note this is the node position with indexing starting at 0. 

    #%% write mesh to file 
    def findIdirichlet(self):
        """Find the best node for the dirichlet node.

        Returns
        -------
        idirchlet : int 
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
        
        idirichlet_edge = np.argmax(np.sqrt(np.sum((edge_nodes - e)**2, axis=1)))
        idirichlet_node = edge_nodes[idirichlet_edge]
        
        idirichlet = np.argmin(np.sqrt(np.sum((self.node - idirichlet_node)**2, axis=1)))

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
            
        # find furthest node from any electrode to be dirichlet node
        self.idirichlet = 0
        if self.idirichlet is None:
            idirichlet = self.findIdirichlet() + 1 
        else: 
            idirichlet = self.idirichlet + 1 
            
        ### write data to mesh.dat kind of file ###
        #open mesh.dat for input      
        with open(file_path, 'w') as fid:
            #write to mesh.dat total num of elements and nodes

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
            
        # element parameters 
        if 'param' in self.df.columns: 
            param = np.asarray(self.df['param'],dtype=int)
        else:
            param = np.arange(self.numel) + 1 
        
        if 'zones' in self.df.columns: 
            zone = np.asarray(self.df['zones'],dtype=self.dint)
        else: 
            zone = np.ones(self.numel,dtype=self.dint) 
        
        nzones = len(np.unique(zone))
            
        adv_flag = int(iadvanced)
        #compute neighbourhood matrix 
        if self.neigh_matrix is None:
            self.computeNeigh()
        neigh = self.neigh_matrix.copy()
        neigh = mc.sortNeigh(neigh,zone) # organise for (c)R3t input  
        neigh += 1 # add one for FOTRAN indexing 
        
        if self.NsizeA is None:#then the finite element conductance matrix needs calculating 
            self.computeNconnec()

        fconm = self.fconm.copy() + 1 #(add 1 for FORTRAN indexing)

        # find furthest node from first electrode to be dirichlet node
        if self.idirichlet is None:
            idirichlet = self.findIdirichlet() + 1 
        else: 
            idirichlet = self.idirichlet + 1 
            
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
                fid.write('{:d}\n'.format(idirichlet))
                
            if nzones > 1: 
                for i in range(nzones):
                    fid.write('{:d} {:f}\n'.format(i+1,1.0))
                    

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
            if key in ['X','Y','Z','elm_id','cellType']: 
                continue 
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
        
    def csv(self,file_name='mesh.csv'):
        """ Write a .csv file of the mesh, the first 3 columns are the element 
        centres at coordinates x,y,z and the rest of the columns are the 
        attributes in the df

        Parameters
        ----------
        file_name : String, optional
            The default is 'mesh.csv'.

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
        
    def xyz(self,file_name='mesh.xyz',coordParam=None):
        """ Write a .xyz file of the mesh, the first 3 columns are the element 
        centres at coordinates x,y,z and the rest of the columns are the 
        attributes in the df

        Parameters
        ----------
        file_name : String, optional
            The default is 'mesh.csv'.

        """
        if isinstance(file_name,str)==False:
            raise NameError("file_name argument must be a string")

        if not file_name.lower().endswith('.xyz'):
            file_name += '.xyz'
            
        # setup coordinate parameters if not given 
        if coordParam is None: 
            coordParam = {}
            coordParam['a'] = 0 
            if self.elec is not None: 
                if self.iremote is not None:
                    iremote = self.iremote 
                else: 
                    iremote = [False]*self.elec.shape[0]
                coordParam['x0'] = np.min(self.elec[:,0][~iremote])
                coordParam['y0'] = np.min(self.elec[:,1][~iremote])
            else:
                coordParam['x0'] = np.min(self.node[:,0])
                coordParam['y0'] = np.min(self.node[:,1])
                
        self.cellCentres() 
        x_coords=self.df.X.values#get element coordinates
        y_coords=self.df.Y.values
        z_coords=self.df.Z.values 
        df = self.df.copy().reset_index()
        keys0= self.df.keys()
        keys=[]
        #ignore the x y z columns if already in df
        ignore = ['X','Y','Z','param','elm_id','cellType','region']
        for k in keys0:
            if k not in ignore:
                keys.append(k)
        
        #open file 
        fh = open(file_name,'w')
        fh.write('X, Y, Z') # write xyz headers 
        [fh.write(', '+key) for key in keys] # write attribute headers 
        fh.write('\n') # drop a line 
        for i in range(self.numel):
            line = '{:f},{:f},{:f}'.format(x_coords[i],y_coords[i],z_coords[i])
            for key in keys:
                line += ',{:}'.format(df[key][i])
            line += '\n'
            fh.write(line)
        # close file     
        fh.close()
        
        # write out voxel information 
        if self.nvoxel['x'] is None: 
            return # exit function if nvoxel not populated 
        
        voxel_info_file = file_name.replace('.xyz','_voxel_info.txt')
        fh = open(voxel_info_file,'w')
        for key in self.nvoxel.keys(): 
            fh.write('%s : %i\n'%(key,int(self.nvoxel[key])))
        fh.close() 
        
        
    def toVoxelMesh(self, elec=None, elec_type=None, **kwargs):
        """
        Returns a special kind of mesh optimised for visualasation inside of 
        GIS software like Geovisionary. 

        Parameters
        ----------
        elec : array like, pd.DataFrame, optional
            Electrode coordinates. The default is None.
        elec_type : list, optional
            List of electrode types, i.e. 'surface' or 'buried'. The default is
            None (all electrodes are surface electrodes).

        Raises
        ------
        Exception
            If no electrode coordinates provided. 

        Returns
        -------
        mesh : class
            Mesh object with voxel elements 

        """
        if elec is None:
            elec = self.elec 
            
        if elec is None:
            raise Exception('need to provide some electrodes to produce a voxel mesh')
        try: 
            elec_x = elec[:,0]
            elec_y = elec[:,1]
            elec_z = elec[:,2]
        except: 
            elec_x = elec.x.values 
            elec_y = elec.y.values 
            elec_z = elec.z.values 
        if elec_type is None:
            elec_type = self.elec_type 
        
        if 'surface_refinement' not in kwargs.keys() and self.ndims==3:
            if self.surfaceMesh is not None: 
                kwargs['surface_refinement'] = self.surfaceMesh.node.copy() 

        mesh = voxelMesh(elec_x, elec_y, elec_z, elec_type, **kwargs, 
                         force_regular=True)
        mesh.cellCentres() 
        self.cellCentres()
        interp_method = 'triangulate'
        if all(self.df.Y == self.df.Y[0]):
            interp_method = 'nearest'
        elif all(self.df.X == self.df.X[0]):
            interp_method = 'nearest'

        for column in self.df.columns:
            if column in ['X','Y','Z','param','elm_id','cellType','region']: 
                continue 
            if interp_method == 'triangulate': 
                inew = interp.triangulate3d(mesh.df.X.values, mesh.df.Y.values, mesh.df.Z.values, 
                                            self.df.X.values, self.df.Y.values, self.df.Z.values, 
                                            self.df[column].values)
            else: 
                inew = interp.nearest3d(mesh.df.X.values, mesh.df.Y.values, mesh.df.Z.values, 
                                        self.df.X.values, self.df.Y.values, self.df.Z.values, 
                                        self.df[column].values)
            mesh.addAttribute(inew, column)
            
        return mesh
    
    
    def vts(self,file_name='mesh.vts',x0=None,y0=None,a=0):
        """
        Structured grid file, therefore must be a Voxel type mesh. Note that 
        national grid parameters are only commented in the resulting vts file 
        for now. 

        Parameters
        ----------
        file_name : str, optional
            Path to output file. The default is 'mesh.vts'.
        x0 : float, optional
            Grid origin in X axis. The default is None.
        y0 : float, optional
            Grid origin in Y axis. The default is None.
        a : float, optional
            Rotation angle or bearing in X direction. The default is 0.

        Returns
        -------
        ~.vts: file 
            vts file written to file_name location. 

        """
        if self.cell_type[0] != 11:
            warnings.warn('.vts export only available for voxel mesh, function call terminated...')
            return 
        
        # write out voxel information 
        if self.nvoxel['x'] is None: 
            warnings.warn('.vts export only available for voxel mesh (with nx, ny, nz information), function call terminated...')
            return # exit function if nvoxel not populated 
        
        fh = open(file_name,'w')
        
        if x0 is not None and y0 is not None: 
            # make a comment line 
            fh.write('<!-- rotation info: x0 = %f; y0 = %f; a = %f -->\n'%(x0,y0,a))
        
        def writeXMLline(text,bracket=True,tab=0):
            for i in range(tab):
                fh.write('\t')
            if bracket: 
                fh.write('<')
            fh.write('{:}'.format(text))
            if bracket:
                fh.write('>')
            fh.write('\n')
            
        writeXMLline('VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian"')
    
        voxel_info = '"0 %i 0 %i 0 %i"'%(self.nvoxel['x'], self.nvoxel['y'], self.nvoxel['z'])

        writeXMLline('StructuredGrid WholeExtent='+voxel_info,tab=1)
        writeXMLline('Piece Extent='+voxel_info, tab=1)
        
        # TODO: write out point data?? 
        writeXMLline('PointData', tab=2)
        writeXMLline('/PointData', tab=2)
        
        # write cell data 
        writeXMLline('CellData',tab=2)
        for column in self.df.columns:
            if column in ['X','Y','Z','param','elm_id','cellType','region']: 
                continue 
            header = 'DataArray type="Float32" Name="%s" format="ascii"'%column 
            writeXMLline(header,tab=3)
            X = self.df[column].values 
            X[np.isnan(X)] = -9999
            for i in range(len(X)):
                writeXMLline('%f'%X[i], False, 4)
            writeXMLline('/DataArray',tab=3)
        writeXMLline('/CellData',tab=2)
        
        # write out point data
        writeXMLline('Points',tab=2)
        writeXMLline('DataArray type="Float32" Name="Points" NumberOfComponents="3" format="ascii"',tab=3)
        for i in range(self.numnp):
            line = '{:f}\t{:f}\t{:f}'.format(self.node[i,0], self.node[i,1], self.node[i,2])
            writeXMLline(line, False, 4)
        writeXMLline('/DataArray',tab=3)
        writeXMLline('/Points',tab=2)
        writeXMLline('/Piece',tab=1)
        writeXMLline('/StructuredGrid',tab=1)
        writeXMLline('/VTKFile')
        
        fh.close()
        
    def vtu(self,file_name='mesh.vtu'):
        """
        Unstructured xml type grid file. 
        
        Parameters
        ----------
        file_name : str, optional
            Path to output file. The default is 'mesh.vtu'.
   
        Returns
        -------
        ~.vts: file 
            vts file written to file_name location. 
   
        """
        fh = open(file_name,'w')
        
        def writeXMLline(text,bracket=True,tab=0):
            for i in range(tab):
                fh.write('\t')
            if bracket: 
                fh.write('<')
            fh.write('{:}'.format(text))
            if bracket:
                fh.write('>')
            fh.write('\n')
            
        def writeXMLarray(X):
            # tab = 5
            j = 0 
            text = '' 
            for i in range(len(X)): 
                if isinstance(X[i],int):
                    text += '{:d} '.format(X[i])
                else: 
                    text += '{:f} '.format(X[i])
                j += 1 
                if j == 5: 
                    writeXMLline(text,False,5)
                    j = 0 
                    text = ''  
            if j < 5: # write out last line 
                writeXMLline(text,False,5)
                
        
        data_header_template = 'DataArray type="Float64" Name="{:s}" format="ascii" RangeMin="{:f}" RangeMax="{:f}"'
        point_header_template = 'DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii" RangeMin="{:f}" RangeMax="{:f}"'
        cell_header_template = 'DataArray type="Int64" Name="connectivity" format="ascii" RangeMin="0" RangeMax="{:f}"'
        ofst_header_template = 'DataArray type="Int64" Name="offsets" format="ascii" RangeMin="{:d}" RangeMax="{:d}"'
        type_header_template = 'DataArray type="UInt8" Name="types" format="ascii" RangeMin="{:d}" RangeMax="{:d}"'
       
        writeXMLline('VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64"')
        writeXMLline('UnstructuredGrid',tab=1)
        writeXMLline('Piece NumberOfPoints="%i" NumberOfCells="%i"'%(self.numnp, self.numel), tab=2)
        
        # write out point data
        writeXMLline('PointData', tab=3)
        for column in self.ptdf.columns:
            X = self.ptdf[column].values 
            header = data_header_template.format(column, np.min(X), np.max(X))
            writeXMLline(header, tab=4)
            X[np.isnan(X)] = -9999
            writeXMLarray(X)
            writeXMLline('/DataArray',tab=4)
        writeXMLline('/PointData', tab=3)
        
        # write cell data 
        writeXMLline('CellData',tab=3)
        for column in self.df.columns:
            if column in ['X','Y','Z','cellType']: 
                continue 
            X = self.df[column].values 
            header = data_header_template.format(column, np.min(X), np.max(X))
            writeXMLline(header,tab=4)
            X[np.isnan(X)] = -9999
            writeXMLarray(X)
            writeXMLline('/DataArray',tab=4)
        writeXMLline('/CellData',tab=3)
        
        # write out node/point coordinates 
        writeXMLline('Points',tab=3)
        header = point_header_template.format(np.min(self.node.flatten()), np.max(self.node.flatten()))
        writeXMLline(header, tab=4)
        for i in range(self.numnp):
            line = '{:f}\t{:f}\t{:f}'.format(self.node[i,0], self.node[i,1], self.node[i,2])
            writeXMLline(line, False, 5)
        writeXMLline('/DataArray',tab=4)
        writeXMLline('/Points',tab=3)
        
        # write out connection matrix 
        writeXMLline('Cells',tab=3)
        header = cell_header_template.format(self.numel-1)
        writeXMLline(header,tab=4)
        for i in range(self.numel):
            text = ''
            for j in range(self.connection.shape[1]):
                text += '{:d}\t'.format(self.connection[i,j])
            writeXMLline(text,False,5)
        writeXMLline('/DataArray',tab=4)
        # write out offsets (take the number of nodes per element, and write out its times table)
        npere = self.type2VertsNo()
        offsets = [npere]*self.numel 
        for i in range(1, self.numel):
            offsets[i] = (i+1)*npere 
        header = ofst_header_template.format(min(offsets),max(offsets))
        writeXMLline(header,tab=4)
        writeXMLarray(offsets)
        writeXMLline('/DataArray',tab=4)
        
        # dont forget the cell types! 
        cell_type = [int(self.cell_type[0])]*self.numel 
        header = type_header_template.format(cell_type[0], cell_type[0])
        writeXMLline(header,tab=4)
        writeXMLarray(cell_type)
        writeXMLline('/DataArray',tab=4)
        writeXMLline('/Cells',tab=3)
        
        # closing lines 
        writeXMLline('/Piece',tab=2)
        writeXMLline('/UnstructuredGrid',tab=1)
        writeXMLline('/VTKFile')
        
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
        ----------
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
            
    
    def exportTetgenMesh(self,prefix='mesh',zone=None, debug=False, mixed=False):
        """ Export a mesh like the tetgen format for input into E4D. 
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
        mixed: bool
            If true attempt to compute mixed boundary condition at mesh edges. 
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
        if not mixed: 
            node_bd[node_bd==2] = 1 
            face_bd[face_bd==2] = 1 
            
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
        """Save mesh into a file. Available formats are .dat, .vtk, .csv, .xyz,
        and .node. .vts is a special case where the mesh should be transformed 
        into a voxel mesh first. 

        Parameters
        ----------
        fname : str, optional 
            Path to output save file. 
        ftype: str, optional
            File extension, if none, will be guessed from file name.
        """
        if not isinstance(fname,str):
            raise TypeError('fname needs to be a string!')
        
        #determine file type 
        atypes = sorted(['dat','node','vtk','vts','vtu','csv', 'xyz'])
        if ftype is None:#guess
            for a in atypes:
                if fname.endswith('.' + a):
                    ftype = a
                    break
        
        if ftype not in atypes:
            raise NameError('Unregocnised mesh file format(%s), avaialable types are %s'%(ftype,str(atypes)))
        
        #add extension if not already there 
        if not fname.endswith('.' + ftype):
            # print('adding file extension')
            fname = fname + '.' + ftype
            
        #call save file function 
        if ftype == 'dat':
            self.dat(fname)
        elif ftype == 'vtk':
            self.vtk(fname)
        elif ftype == 'vtu':
            self.vtu(fname)
        elif ftype == 'csv':
            self.csv(fname)
        elif ftype == 'node':
            self.exportTetgenMesh(fname.replace('.node',''))
        elif ftype == 'xyz':
            self.xyz(fname) 
        elif ftype == 'vts':
            self.vts(fname)
        else:
            raise ValueError('mesh export format not recognized. Try either .vtk, .node or .dat.')
            
    def exportMesh(self, fname, ftype=None, coordLocal=False, coordParam=None, 
                   voxel=False, meshParam=None):
        """
        Export mesh into other formats, allowing for coordinate rotation.  
        
        Parameters
        ----------
        fname : str, optional 
            Path to output save file. 
        ftype: str, optional
            File extension, if none, will be guessed from file name.
        coordLocal: bool, optional 
            If True, convert mesh coordinates to a national grid system or utm. 
        coordParam: dict, optional
            Coordinate conversion parameters, x0, y0 and a. Stored as a dictionary. 
        """
        # check if vts format! This will force results to be as a voxel and in local grid 
        vtsflag = False 
        if ftype is not None and ftype.endswith('vts'): 
            vtsflag = True 
        if fname.endswith('vts'):
            vtsflag = True 
        if vtsflag: 
            # can't have global coordinates with a vts file 
            if coordLocal: 
                warnings.warn('VTS export will be written in a local grid format!')
            coordLocal = False 
            voxel = True # also vts must be a voxel type format! 
            
        if coordLocal and coordParam is None: 
            coordParam = {}
            coordParam['a'] = 0 
            if self.elec is not None: 
                if self.iremote is not None:
                    iremote = self.iremote 
                else: 
                    iremote = [False]*self.elec.shape[0]
                coordParam['x0'] = np.min(self.elec[:,0][~iremote])
                coordParam['y0'] = np.min(self.elec[:,1][~iremote])
            else:
                coordParam['x0'] = np.min(self.node[:,0])
                coordParam['y0'] = np.min(self.node[:,1])
        
        kwargs = {} 
        kwargs['cl_factor'] = 1 # force mesh to the extent of the electrodes 
        if meshParam is not None: 
            if 'cl' in meshParam.keys():
                kwargs['cl'] = meshParam['cl']
            if 'fmd' in meshParam.keys():
                kwargs['fmd'] = meshParam['fmd']
                
        if self.elec is not None:
            dist = np.unique(findDist(self.elec[:,0], self.elec[:,1], self.elec[:,2]))[1]
            kwargs['cl'] = dist/2 
                
        if voxel: 
            mesh = self.toVoxelMesh(**kwargs)
        else:
            mesh = self.copy() 
            
        if coordLocal: 
            x0 = coordParam['x0']
            y0 = coordParam['y0']
            a = coordParam['a']
            utmx, utmy = interp.invRotGridData(
                mesh.node[:,0], 
                mesh.node[:,1], 
                x0, y0, a)
            mesh.node[:,0] = utmx 
            mesh.node[:,1] = utmy 
            mesh.saveMesh(fname,ftype)
        else:
            mesh.saveMesh(fname,ftype)
        
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
        px = np.polyfit(xLocal, xGrid,1)
        py = np.polyfit(xGrid, yGrid,1)
        mesh.elec[:,0][~mesh.iremote] = np.polyval(px, xLocal) 
        mesh.elec[:,1][~mesh.iremote] = np.polyval(py, xGrid) 
        mesh.node[:,0] = np.polyval(px, mesh.node[:,0])
        mesh.node[:,1] = np.polyval(py, mesh.node[:,0]) 

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
    # read in header lines 
    title=fid.readline().strip()#read line 2
    while title == '':
        title=fid.readline().strip()#read line 2
    
    format_type=fid.readline().strip()#read line 3
    while format_type == '':
        format_type=fid.readline().strip()#read line 3
        
    if format_type=='BINARY':
        raise ImportError("expected ASCII type file format, not binary")
        
    dataset_type=fid.readline().strip().split()#read line 4
    while len(dataset_type) == 0: 
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
    rowwise = True # flag that nodes are listed rowwise 
    for i in range(numnp):
        node_num[i] = i
        coord_data_line = fid.readline().strip()
        coord_data = coord_data_line.split() 
        if len(coord_data) == 3: 
            node_x[i] = float(coord_data[0])
            node_y[i] = float(coord_data[1])
            node_z[i] = float(coord_data[2])
        elif len(coord_data) == 1:
            coord_data=fid.readline()
            node_x[i] = float(coord_data_line[0:12]) # retrive fixed width columns if cannot parse as split strings
            node_y[i] = float(coord_data_line[12:24])
            node_z[i] = float(coord_data_line[24:36])
        else: 
            rowwise = False # switch off flag that nodes are row wise 
            break 
    
    if not rowwise: 
        #read in node data that is in matrix form 
        node = np.zeros((numnp,3),dtype=float)
        i = 0 
        j = 0 
        k = 0 
        while i < len(coord_data):
            node[k,j] = float(coord_data[i])
            j += 1 
            i += 1 
            if j == 3:
                j = 0 
                k += 1 
                
        while k < numnp: 
            coord_data_line = fid.readline().strip()
            coord_data = coord_data_line.split() 
            i = 0 
            j = 0 
            while i < len(coord_data):
                node[k,j] = float(coord_data[i])
                j += 1 
                i += 1 
                if j == 3:
                    j = 0 
                    k += 1 
                
        node_x = node[:,0]
        node_y = node[:,1]
        node_z = node[:,2]
    
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


def vtuImport(file_path,order_nodes=True):
    """Vtk importer for newer format (not fully validated)
    """
    if os.path.getsize(file_path)==0: # So that people dont ask me why you cant read in an empty file, throw up this error. 
        raise ImportError("Provided mesh file is empty! Check that (c)R2/3t code has run correctly!")
        
    # get vtu file information (not robust yet)
    root = ET.parse(file_path)
    for child in root.iter('VTKFile'):
        byte_order = child.attrib['byte_order']
        vtk_type = child.attrib['type']
    
    if vtk_type != 'UnstructuredGrid':
        raise Exception('vtu file of unrecognised grid type, only unstructed grids allowed')

    # function to convert acsii data contained in an xml element array to python list 
    def child2Array(child):
        if 'format' not in child.attrib.keys():
            raise Exception('Cannot find format identifier in vtu file')
        fmt = child.attrib['format']
        if fmt not in ['ascii','binary']: 
            raise Exception('Vtu data array not of ASCII  or Binary type, aborting...')
            
        typ = float
        bcode = '<d' # byte decoding code 
        if byte_order == 'BigEndian':
            bcode = '>d'
        pad = 1 # padding 
        if 'int' in child.attrib['type'].lower():
            typ = int
            bcode = bcode.replace('b','q') # change code to long int 
            if 'uint8' in child.attrib['type'].lower():
                bcode = bcode.replace('q','B') # change code to unsigned char 
                pad = 8 
            
        t = child.text.strip() # get raw text 
        if fmt == 'ascii': 
            t = t.replace('\n','')
            a = [typ(x) for x in t.split()]
        else:
            b64 = bytes(t, encoding='ascii') # convert to ascii with base64 binary 
            b = base64.decodebytes(b64)
            a = [u[0] for u in struct.iter_unpack(bcode,b)]
            a = a[pad:]
        return a 


    # get the number of cells and points 
    for child in root.iter('Piece'):
        numnp = int(child.attrib['NumberOfPoints'])
        numel = int(child.attrib['NumberOfCells'])
        
    # get point data 
    ptdf = pd.DataFrame()
    for child in root.iter('PointData'):
        for subchild in child.findall('DataArray'):
            a = child2Array(subchild)
            ptdf[subchild.attrib['Name']]=a
            
    # get cell data 
    df = pd.DataFrame()
    for child in root.iter('CellData'):
        for subchild in child.findall('DataArray'):
            a = child2Array(subchild)
            df[subchild.attrib['Name']]=a
            
    # get node coordinates 
    for child in root.iter('Points'):
        for subchild in child.findall('DataArray'):
            if not subchild.attrib['Name'] == 'Points':
                continue 
            ncomp = int(subchild.attrib['NumberOfComponents'])
            assert ncomp == 3 
            shape = (numnp,ncomp)
            node = np.array(child2Array(subchild)).reshape(shape)
            
    # get connectivity matrix properties 
    offset = None 
    cshape = None 
    for child in root.iter('Cells'):
        for subchild in child.findall('DataArray'):
            if not subchild.attrib['Name'] == 'offsets':
                continue 
            offset = int(subchild.attrib['RangeMin'])
            cshape = (numel,offset)
            
    if offset is None: 
        raise Exception('Cannot find offsets element in vtu file, cannot build connection matrix! Aborting...')
            
    cellType = None 
    for child in root.iter('Cells'):
        for subchild in child.findall('DataArray'):
            if not subchild.attrib['Name'] == 'types':
                continue 
            cellType = child2Array(subchild)
            
    # now find connectivity matrix 
    for child in root.iter('Cells'):
        for subchild in child.findall('DataArray'):
            if not subchild.attrib['Name'] == 'connectivity':
                continue 
            connection = np.array(child2Array(subchild)).reshape(cshape)
                
    mesh = Mesh(node[:,0],#x coordinates of nodes 
                node[:,1],#y coordinates of nodes
                node[:,2],#z coordinates of nodes 
                connection,#nodes of element vertices
                cellType,#according to vtk format
                file_path,
                order_nodes) 
    
    #add attributes / cell parameters 
    mesh.df = df 
    mesh.ptdf = ptdf 

    return mesh


def vtsImport(file_path,order_nodes=True):
    """Vtk importer for newer format (not fully validated)
    """
    if os.path.getsize(file_path)==0: # So that people dont ask me why you cant read in an empty file, throw up this error. 
        raise ImportError("Provided mesh file is empty!")
        
    # get vtu file information (not robust yet)
    root = ET.parse(file_path)
    for child in root.iter('VTKFile'):
        byte_order = child.attrib['byte_order']
        vtk_type = child.attrib['type']
    
    if vtk_type != 'StructuredGrid':
        raise Exception('vts file of unrecognised grid type, only structed grids allowed')

    # function to convert acsii data contained in an xml element array to python list 
    def child2Array(child):
        if 'format' not in child.attrib.keys():
            raise Exception('Cannot find format identifier in vtu file')
        fmt = child.attrib['format']
        if fmt not in ['ascii','binary']: 
            raise Exception('Vtu data array not of ASCII  or Binary type, aborting...')
            
        typ = float
        bcode = '<d' # byte decoding code 
        if byte_order == 'BigEndian':
            bcode = '>d'
        pad = 1 # padding 
        if 'int' in child.attrib['type'].lower():
            typ = int
            bcode = bcode.replace('b','q') # change code to long int 
            if 'uint8' in child.attrib['type'].lower():
                bcode = bcode.replace('q','B') # change code to unsigned char 
                pad = 8 
            
        t = child.text.strip() # get raw text 
        if fmt == 'ascii': 
            t = t.replace('\n','')
            a = [typ(x) for x in t.split()]
        else:
            b64 = bytes(t, encoding='ascii') # convert to ascii with base64 binary 
            b = base64.decodebytes(b64)
            a = [u[0] for u in struct.iter_unpack(bcode,b)]
            a = a[pad:]
        return a 

    cellType = [12]

    # get the number of cells and points 
    for child in root.iter('Piece'):
        extents = child.attrib['Extent'].split()
        numnpx = int(extents[1])+1
        numnpy = int(extents[3])+1
        numnpz = int(extents[5])+1
        numnp = numnpx * numnpy * numnpz
        numel = int(extents[1])*int(extents[3])*int(extents[5])
        
    # get point data 
    ptdf = pd.DataFrame()
    for child in root.iter('PointData'):
        for subchild in child.findall('DataArray'):
            a = child2Array(subchild)
            ptdf[subchild.attrib['Name']]=a
            
    # get cell data 
    df = pd.DataFrame()
    for child in root.iter('CellData'):
        for subchild in child.findall('DataArray'):
            a = child2Array(subchild)
            df[subchild.attrib['Name']]=a
            
    # get node coordinates 
    for child in root.iter('Points'):
        for subchild in child.findall('DataArray'):
            if not subchild.attrib['Name'] == 'Points':
                continue 
            ncomp = int(subchild.attrib['NumberOfComponents'])
            assert ncomp == 3 
            shape = (numnp,ncomp)
            node = np.array(child2Array(subchild)).reshape(shape)
             
    # get connectivity matrix properties 
    connection = np.zeros((numel, 8), dtype=int) 
    xcounter = 1
    ycounter = 1
    zcounter = 1
    i = 0 
    j = 0 
    while i < numnp: 
        #note: the nodes are sorted by x, y, then z. 
        if xcounter == numnpx: 
            ycounter += 1
            xcounter = 1
            i += 1 # jump to next x row 
        
        if ycounter == numnpy:
            zcounter += 1 
            ycounter = 1
            xcounter = 1
            i += numnpx # jump to next y row 
            
        if zcounter == numnpz:
            # reached top or bottom of mesh, so break here 
            break 
        
        connection[j, 0] = i
        connection[j, 1] = i + 1
        connection[j, 4] = i + numnpx 
        connection[j, 5] = i + numnpx + 1 
        connection[j, 3] = i + (numnpy * numnpx)
        connection[j, 2] = i + (numnpy * numnpx) + 1 
        connection[j, 7] = i + (numnpy * numnpx) + numnpx
        connection[j, 6] = i + (numnpy * numnpx) + numnpx + 1  
        j += 1 
        i += 1 
        xcounter += 1 
                
    mesh = Mesh(node[:,0],#x coordinates of nodes 
                node[:,1],#y coordinates of nodes
                node[:,2],#z coordinates of nodes 
                connection,#nodes of element vertices
                cellType,#according to vtk format
                file_path,
                False) 
    
    #add attributes / cell parameters 
    mesh.df = df 
    mesh.ptdf = ptdf 

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
        zone[i] = int(npere+2)
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
    order_nodes: bool, optional
        Check ordering of the nodes on import, this ensures that the mesh will 
        work inside of R3t / cR3t. Default is True. 
    
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
    
    mesh.addAttribute(zone,'zones')
    
    return mesh
        
        
#%% build a quad mesh        
def quadMesh(elec_x, elec_z, elec_type = None, elemx=4, xgf=1.5, zf=1.1, zgf=1.5, 
             fmd=None, pad=2, surface_x=None, surface_z=None,
             refine_x = None, refine_z=None, model_err=False):
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
    model_err : bool, optional 
        If True all topography will be normalised to zero such that the returned 
        mesh can be used to model foward modelling errors. 
        
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
    
    elec_x = np.array(elec_x)
    elec_z = np.array(elec_z)
    elec_x_cache = None 
    elec_z_cache = None 
    
    bh_flag = False
    surface_idx = []#surface electrode index
    bur_idx = []#buried electrode index 
    rem_idx = [] 
    #determine the relevant node ordering for the surface electrodes? 
    if elec_type is not None:
        if not isinstance(elec_type,list):
            raise TypeError("electrode_type variable is given but expected a list type argument")
        if len(elec_type) != len(elec_x):
            raise ValueError("mis-match in length between the electrode type and number of electrodes")
        
        for i,key in enumerate(elec_type):
            if key == 'electrode': 
                surface_idx.append(i)
            if key == 'buried': 
                bur_idx.append(i)
                bh_flag = True
            if key == 'remote':
                rem_idx.append(i)
            
        if len(surface_idx)>0:# then surface electrodes are present
            Ex = np.array(elec_x)[surface_idx]
            Ez = np.array(elec_z)[surface_idx]
        elif len(surface_idx) == 0 and len(surface_x) > 0:
            #case where you have surface topography but no surface electrodes 
            Ex = np.array(surface_x)
            Ez = np.array(surface_z)
            elec = np.c_[Ex, Ez]
        elif len(surface_idx)== 0:
            #fail safe if no surface electrodes are present to generate surface topography 
            Ex = np.array([elec_x[np.argmin(elec_x)], elec_x[np.argmax(elec_x)]])
            Ez = np.array([elec_z[np.argmax(elec_z)], elec_z[np.argmax(elec_z)]])
    else:
        surface_idx = np.arange(len(elec_x))
        
    elec = np.c_[elec_x,elec_z]    
    
    if bh_flag:
        bh = np.c_[np.array(elec_x)[bur_idx],np.array(elec_z)[bur_idx]]
        
    if model_err: 
        # if attempting to model forward modelling errors, normalise topography
        x_interp = np.append(elec[:,0][surface_idx],surface_x)
        z_interp = np.append(elec[:,1][surface_idx],surface_z)
        ifunc = interp1d(x_interp, z_interp, fill_value='extrapolate')
        elec_z = np.array(elec_z) - ifunc(np.array(elec_x))
        elec = np.c_[elec_x,elec_z]  
        if len(surface_z) > 0:
            surface_z = surface_z - ifunc(surface_x)
        if elec_type is not None: 
            Ez = Ez - ifunc(Ex)
        if refine_flag: 
            refine_z = np.array(refine_z) - ifunc(np.array(refine_x))
        if bh_flag:
            bh[:,1] = bh[:,1] - ifunc(bh[:,0])
            
    if len(rem_idx) > 0:
        elec_x_cache = elec_x.copy()
        elec_z_cache = elec_z.copy()
        elec_x = np.delete(elec_x,rem_idx)
        elec_z = np.delete(elec_z,rem_idx)
        elec = np.c_[elec_x,elec_z]
        
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
        
    meshz = [0]
    dyy = np.mean(np.diff(elecXsorted))/5
    for i in range(100):
        meshz.append(meshz[-1]+dyy*zf)
        dyy = dyy*zf
        if meshz[-1] > fmd:
            break
    elemy = len(meshz)
    elemy2 = int(elemy*1.2)
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
        Y = np.append(Ez, surface_z)
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
    
    ## find the columns relating to the electrode nodes? ##
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
    if len(rem_idx) > 0:
        # get original x z coordinates with the remote electrode included 
        elec_x = elec_x_cache.copy() 
        elec_z = elec_z_cache.copy() 
        for n,i in enumerate(rem_idx):
            elec_x[i] = meshx[n]
            elec_z[i] = np.min(node_z)
        
    node_in_mesh = [0]*len(elec_x)
    for i in range(len(elec_x)):
        sq_dist = (node_x - elec_x[i])**2 + (node_z - elec_z[i])**2 # find the minimum square distance
        
        node_in_mesh[i] = np.argmin(sq_dist) # min distance should be zero, ie. the node index.
    mesh.setElecNode(node_in_mesh) # add nodes to the mesh class
    mesh.elec_type = elec_type
    
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
        
    if dump is print: 
        def dump(x):
            print(x,end='')

    if platform.system() == "Windows":#command line input will vary slighty by system 
        cmd_line = [os.path.join(ewd,'gmsh.exe'), file_name+'.geo', opt, 'nt %i'%ncores]
        
    elif platform.system() == 'Darwin': # its a macOS 
        # winetxt = 'wine'
        # if getMacOSVersion():
        #     winetxt = 'wine64'
        winetxt = whichWineMac()
        winePath = []
        wine_path = Popen(['which', winetxt], stdout=PIPE, shell=False, universal_newlines=True)#.communicate()[0]
        for stdout_line in iter(wine_path.stdout.readline, ''):
            winePath.append(stdout_line)
        if winePath != []:
            cmd_line = ['%s' % (winePath[0].strip('\n')), ewd+'/gmsh.exe', file_name+'.geo', opt,'-nt','%i'%ncores]
        else:
            cmd_line = [wPath, ewd+'/gmsh.exe', file_name+'.geo', opt,'-nt','%i'%ncores]   
        
    elif platform.system() == 'Linux':
        #if platform.machine() == 'aarch64':
        #    cmd_line = [ewd + '/gmsh_aarch64', file_name + '.geo', opt,'-nt','%i'%ncores]
        if platform.machine() in ['armv7l', 'aarch64']:
            cmd_line = [ewd + '/gmsh_armv7l', file_name + '.geo', opt,'-nt','%i'%ncores]
        elif os.path.isfile(os.path.join(ewd,'gmsh_linux')):
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
            line = p.stdout.readline()
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
    if any(dist[:,1]<1e-16): # raise an issue when repeated nodes are present 
        pblm_idx = np.arange(pts.shape[0])[dist[:,1] < 1e-16]
        
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
        error = 'The following nodes are repeated!\n'
        count = 0 
        for pair in pairs:
            for p in pair: 
                x_ = X[p]
                y_ = Y[p]
                z_ = Z[p]
                f_ = flag[p]
                line = "X = {:16.8f}, Y = {:16.8f}, Z = {:16.8f}, flag = {:}\n".format(x_,y_,z_,f_)
                # print(line.strip())
                error += line 
            count += 1 
            if count > 10:
                break # otherwise error message will be too big 
                
        
        #warnings.warn(error,Warning)
        raise ValueError(error)

#%% control point helper functions 

def bearing(dx,dy): # compute bearing 
    if dx == 0 and dy == 0:
        raise ValueError('both dx and dy equal 0 - check no 2 electrodes occupy same xy coordinates')
    elif dx == 0 and dy > 0:
        return 0
    elif dx == 0 and dy < 0:
        return 180
    elif dx > 0 and dy == 0:
        return 90
    elif dx < 0 and dy == 0:
        return 270
    elif dx > 0 and dy > 0: 
        return np.rad2deg(np.arctan(dx/dy))
    elif dx > 0 and dy < 0: 
        return 180 + np.rad2deg(np.arctan(dx/dy))
    elif dx < 0 and dy < 0: 
        return 180 + np.rad2deg(np.arctan(dx/dy))
    elif dx < 0 and dy > 0: 
        return 360 + np.rad2deg(np.arctan(dx/dy))
    
def quad(dx,dy):
    b = bearing(dx,dy)
    if b > 315 or b <= 45:
        return 0 
    elif b > 45 and b <= 135: 
        return 1 
    elif b > 135 and b <= 225: 
        return 2 
    else:
        return 3 
    
def halfspaceControlPoints(elec_x, elec_y, r, min_r=None, check_quadrant=True, 
                           cfactor=5, nfactor=100):
    """
    Create control points for surface electrodes in 3D (or 2d borehole electrodes), 
    which attempts to better refine the mesh further away from the electrodes. 
    
    Parameters
    ----------
    elec_x: array like 
        Electrode X coordinates 
    elec_y: array like 
        Electrode Y coordinates (or Z coordinates in the case of borehole surveys)
    r: float
        Radius of control points. Sets the minimum distance from the electrodes. 
    min_r: float, optional 
        Radius in which multiple control points cannot exist, prevents over tuning
        of mesh refinement near the electrodes. Set at the desired characteristic
        length. Defaults to the same as 'r'. 
    check_quadrant: bool, optional 
        Checks if refinement points are being added inline with electrodes. For
        most cases best to leave this as True. 
    cfactor: float, int, optional
        Characteristic length multiplication factor for near electrode field. 
    nfactor: float, int, optional 
        Neuman region characteristic length factor. Helps enlarge points inserted 
        via Voroni calculations further away from the electrodes. Default is 10. 
    
    Returns:
    --------
    xstack: array like 
        Array of refinement point X coordinates 
    ystack: array like 
        Array of refinement point Y coordinates 
    ptsize: array like 
        Array of suggested characteristic lengths for the refinement points 
    """
    if min_r is None:
        min_r = r
        
    # electrode geometry parameters 
    elec_x = np.array(elec_x)
    elec_y = np.array(elec_y)
    nelec = len(elec_x)
    points = np.c_[elec_x,elec_y]
    tree = cKDTree(points)
    _,idx = tree.query(points,3)
    
    xstack = [] 
    ystack = [] 
    ptsize = [] 
    
    addry = [r,0,-r,0]
    addrx = [0,r,0,-r]
    
    if nelec < 4:
        check_quadrant = False 
    
    if check_quadrant: 
        for i in range(nelec):
            dx0 = elec_x[idx[i,1]] - elec_x[i] 
            dx1 = elec_x[idx[i,2]] - elec_x[i]
            dy0 = elec_y[idx[i,1]] - elec_y[i]
            dy1 = elec_y[idx[i,2]] - elec_y[i]
            
            q0 = quad(dx0, dy0) # quadrant 
            q1 = quad(dx1, dy1) # quadrant 
            qs = [q0,q1]
            
            # avoid adding points in direction of electrode array / line 
            # priorities offline control points 
            for j in range(4):
                if j not in qs: 
                    xstack.append(elec_x[i] + addrx[j]) 
                    ystack.append(elec_y[i] + addry[j]) 
                    ptsize.append(min_r*cfactor)
                    
        # do some veroni points further away from electrodes 
        # but only keep those points which are further away in between electrode lines 
        vor = Voronoi(points)
        verts = vor.vertices.copy() 
        # check points are not close too electrodes 
        closecheck = cKDTree(points).query_ball_point(verts, r) 
        # check points are not too close to each other 
        repetcheck = cKDTree(np.vstack([points,verts])).query_ball_point(verts, r*3) 
        # check points are inside convex hull 
        chull = ConvexHull(points)
        chullpath = mpath.Path(points[chull.vertices])
        
        #loop through to satisfy above checks 
        for i in range(verts.shape[0]):
            keep = True  
            if len(closecheck[i]) > 1: 
                c = 0 
                for c,j in enumerate(closecheck[i]):
                    if c == 0: 
                        continue 
                    keep = False 
            if len(repetcheck[i]) > 1: 
                c = 0 
                for c,j in enumerate(repetcheck[i]):
                    if c == 0: 
                        continue 
                    keep = False  
            if not chullpath.contains_point(verts[i,:]):
                keep = False  
            if keep: 
                xstack.append(verts[i,0])
                ystack.append(verts[i,1])
                ptsize.append(min_r*nfactor)
    else: 
        for i in range(nelec):
            for j in range(4):
                xstack.append(elec_x[i] + addrx[j]) 
                ystack.append(elec_y[i] + addry[j]) 
                ptsize.append(min_r*cfactor)
    
    xstack = np.array(xstack)
    ystack = np.array(ystack)
    ptsize = np.array(ptsize)
    # double check no points are repeated 
    tree = cKDTree(np.c_[xstack, ystack])
    ilookup = tree.query_ball_point(np.c_[xstack, ystack],min_r)
    
    # keep only first instance of points close to each other 
    repeats = np.array([False]*len(xstack),dtype=bool)
    for i in range(xstack.shape[0]):
        if len(ilookup[i]) > 1: 
            c = 0 
            for c,j in enumerate(ilookup[i]):
                if c == 0: 
                    continue 
                repeats[j] = True 
    return xstack[~repeats], ystack[~repeats], ptsize[~repeats] 

#%% build a triangle mesh - using the gmsh wrapper
def triMesh(elec_x, elec_z, elec_type=None, geom_input=None, keep_files=True, 
             show_output=True, path='exe', dump=print, whole_space=False, 
             model_err = False, handle=None, **kwargs):
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
        Flag for if the problem should be treated as a whole space porblem, in which case 
        electrode type is ingored and all electrodes are buried in the middle of a large mesh. 
    model_err : bool, optional 
        If True all topography will be normalised to zero such that the returned 
        mesh can be used to model foward modelling errors. 
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
        
    if 'cl' in kwargs.keys():
        cl = kwargs['cl']
    else:
        dist_sort = np.unique(gw.find_dist(elec_x,np.zeros_like(elec_x),elec_z))
        cl = dist_sort[1]/2 # characteristic length is 1/2 the minimum electrode distance
        
    if 'cl_factor' in kwargs.keys():
        cl_factor = kwargs['cl_factor']
    else: 
        cl_factor = 2.0 

    bur_idx=[]#buried electrode index 
    surface_idx = [] 
    bu_flag = False 
    if elec_type is not None: 
        for i,key in enumerate(elec_type):
            if key == 'buried': 
                bur_idx.append(i); bu_flag = True 
            elif key == 'electrode' or key == 'surface':
                surface_idx.append(i)
    else:
        surface_idx = np.arange(len(elec_x))
        
    if bu_flag: 
        bux = np.array(elec_x)[np.array(bur_idx)]
        buz = np.array(elec_z)[np.array(bur_idx)]
        # if all buried electrodes are on a line, that won't work
        if len(np.unique(bux)) > 2 and len(np.unique(buz)) > 2:
            # control points in x z coordinates 
            cpx, cpz, cpl = halfspaceControlPoints(bux, buz, 
                                                   cl, cl*1.5, cl_factor)
            if geom_input is None: 
                geom_input = {'refine':[cpx,cpz,cpl]} 
            elif 'refine' in geom_input.keys():
                geom_input['refine'][0] = list(geom_input['refine'][0]) + cpx.tolist()
                geom_input['refine'][1] = list(geom_input['refine'][1]) + cpz.tolist()
            else:
                geom_input['refine'] = [cpx,cpz,cpl]        
            
    if model_err and not whole_space: 
        # if attempting to model forward modelling errors, normalise topography
        x_interp = np.array(elec_x)[surface_idx]
        z_interp = np.array(elec_z)[surface_idx]
        if geom_input is not None: 
            if 'surface' in geom_input.keys(): 
                x_interp = np.append(x_interp, geom_input['surface'][0])
                z_interp = np.append(z_interp, geom_input['surface'][1])
        ifunc = interp1d(x_interp, z_interp, fill_value='extrapolate')
        elec_z = np.array(elec_z) - ifunc(np.array(elec_x))
        if geom_input is not None: 
            for key in geom_input.keys(): 
                geom_input[key][1] = np.array(geom_input[key][1]) - ifunc(np.array(geom_input[key][0]))
    
    #make .geo file
    file_name="mesh"
    if not whole_space:#by default create survey with topography 
        node_pos = gw.halfspace2d([elec_x,elec_z], elec_type, geom_input,
                                  file_path=file_name,**kwargs)
    elif whole_space:
        print("Whole space problem")
        node_pos = gw.wholespace2d([elec_x,elec_z], 20, elec_type, 
                                   geom_input = geom_input, 
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
    mesh.elec_type = elec_type
    
    return mesh

#%% 3D tetrahedral mesh 
def tetraMesh(elec_x, elec_y, elec_z=None, elec_type = None, keep_files=True, 
              interp_method = 'triangulate', surface_refinement=None, 
              mesh_refinement=None, ball_refinement=True,
              path='exe', dump=print, whole_space=False, model_err=False,
              handle=None, show_output=True, **kwargs):
    """Generates a tetrahedral mesh for R3t (with topography) for field surveys.
    This function expects the current working directory has path: exe/gmsh.exe.
    Uses post processing after mesh generation to super impose topography on to 
    a flat 3D tetrahedral mesh of a hemisphere. Handles, wholespace, 2d and 
    halfspace problems. 
            
    Parameters
    ---------- 
    elec_x : array like
        electrode x coordinates 
    elec_y : array like 
        electrode y coordinates 
    elec_z : array like 
        electrode z coordinates 
    elec_type : list of strings, optional
        Defines if electrodes are buried or not.   
    keep_files : boolean, optional
        `True` if the gmsh input and output file is to be stored in the working 
        directory.
    interp_method : string, default ='triangulate' optional
        Interpolation method to translate mesh nodes in the z direction. In 
        other words the method in which topography is appended to the mesh. 
        Here the topography is added to the mesh in post processing. 
        The options are documented below in the notes. 
    surface_refinement : np.array, optional 
        Numpy array of shape (3,n), should follow the format 
        np.array([x1,x2,x3,...],[y1,y2,y3,...],[z1,z2,z3,...]).
        Allows for extra refinement for the top surface of the mesh. 
        The points are not added to the mesh, but considered in post processing 
        in order to super impose topography on the mesh. 
    mesh_refinement : dict, pd.DataFrame, optional 
        Dataframe (or dict) contianing 'x', 'y', 'z', 'type', and 'cl', columns 
        which describe control points which can be used to refine the mesh, unlike 
        surface_refinement, this argument allows the user granular control over 
        the refinement of mesh elements. If not provided ResIPy attempts to add 
        mesh_refinement for you (though not in the case of a wholespace). 
        See further explanation in tetraMesh notes. 
    ball_refinement : boolean, optional
        If True, tells gmsh to add a 'ball' of refined mesh around electrodes. 
        Default is True. 
    path : string, optional
        Path to exe folder (leave default unless you know what you are doing).
    whole_space : boolean, optional
        flag for if the problem should be treated as a whole space porblem, in 
        which case electrode type is ignored and all electrodes are buried in 
        the middle of a large mesh.
    model_err : bool
        If True, a flat mesh will be returned for the sake of estimating 
        forward modelling errors. 
    dump : function, optional
        Function to which pass the output during mesh generation. `print()` is
        the default.
    **kwargs : dict
        Keyword arguments to be passed to functions in gmshWrap.py 
         
    Returns
    -------
    mesh3d: class
    
    
    Notes 
    -----
    Possible arguments for interp_method: 
        'bilinear' : 4 known points are used to compute the equation of a plane 
                    in which the interpolated point lies. This method is reccomended 
                    if elevation data is organised in a regular grid. 
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
                    
    Format of mesh_refinement:
        The variable must have x,y,z and type columns
        'x': x coordinate array like 
        'y': y coordinate array like 
        'z': z coordinate array like 
        'type': list, object array of point types, use the tag 'surface' for 
            surface points, use 'buried' for buried points. 
        'cl': array like of characteristic lengths for each refinement point in 
            the mesh. 
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
    avail_methods = ['bilinear','nearest','spline','triangulate',None]
    if interp_method not in avail_methods:
        raise NameError("'%s' is an unrecognised interpretation method"%interp_method)
    if elec_z is None: 
        elec_z = np.zeros_like(elec_x)
    if all(np.array(elec_z)==0):
        interp_method = None 
        
    lineis2d = False 
    if all(np.array(elec_y)==elec_y[0]) and 'buried' not in elec_type:
        # do 2d line in 3D mesh 
        lineis2d = True 
        if len(elec_x) > 40: #somewhat arbitary threshold ??
            # then its a long line , use a cylinderical mesh 
            if 'geom' not in kwargs.keys():
                kwargs['geom'] = 'cy'
            if 'flank_fac' not in kwargs.keys():
                kwargs['flank_fac'] = 10 
                
    coplaner = False  
    if all(np.array(elec_y)==elec_y[0]) and surface_refinement is None:
        # cant use triangulation methods if there is no surface refinement and 
        # there is no depth in the Y axis 
        interp_method = 'nearest'
        coplaner = True 
    elif all(np.array(elec_x)==elec_x[0]) and surface_refinement is None:
        # cant use triangulation methods if there is no surface refinement and 
        # there is no depth in the X axis 
        interp_method = 'nearest'
        coplaner = True 
        
    if whole_space: # if whole space problem ignore interpolation 
        interp_method = None 
        
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
    surf_elec_x = []
    surf_elec_y = []
    surf_elec_z = []
    bur_elec_x = []
    bur_elec_y = []
    bur_elec_z = []
    surf_elec_idx = []
    bur_elec_idx = []
    if elec_type is not None and whole_space==False:
        if not isinstance(elec_type,list):
            raise TypeError("'elec_type' argument should be of type 'list', got type %s"%str(type(elec_type)))
        elif len(elec_type) != len(elec_x):
            raise ValueError("mismatch in elec_type vector and elec_x vector lengths")
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
        elec_z = np.array(elec_z) 
        
        if interp_method == 'triangulate':
            # need to check the number of interpolation points is stable for triangulation 
            if len(x_interp) == 4: 
                interp_method = 'bilinear'
            elif len(x_interp) < 4: 
                interp_method = 'nearest'
        
        if len(bur_elec_x)>0: #if we have buried electrodes normalise their elevation to as if they are on a flat surface
            dump('found buried electrodes')
            if interp_method is None:
                bur_elec_z_topo= np.zeros_like(bur_elec_idx) # still need to normalise electrode depths if we want a flat mesh, so use biliner interpolation instead
            elif interp_method == 'bilinear': 
                bur_elec_z_topo = interp.interp2d(bur_elec_x, bur_elec_y, x_interp, y_interp, z_interp)
            elif interp_method == 'nearest':
                bur_elec_z_topo = interp.nearest(bur_elec_x, bur_elec_y, x_interp, y_interp, z_interp)
            elif interp_method == 'spline':
                bur_elec_z_topo = interp.interp2d(bur_elec_x, bur_elec_y, x_interp, y_interp, z_interp,method='spline')
            elif interp_method == 'triangulate':
                bur_elec_z_topo = interp.triangulate(bur_elec_x, bur_elec_y, x_interp, y_interp, z_interp)
    
            elec_z[bur_elec_idx] = elec_z[bur_elec_idx] - bur_elec_z_topo # normalise to zero surface 
            bur_elec_z = elec_z[bur_elec_idx] 

        elec_z[surf_elec_idx] = 0
        
    elif elec_type is not None and whole_space==True:
        if elec_z is None:
            elec_z = np.zeros_like(elec_x)
        for i, key in enumerate(elec_type):
            bur_elec_x.append(elec_x[i])
            bur_elec_y.append(elec_y[i])
            bur_elec_z.append(elec_z[i])
        
    elif not whole_space:
        surf_elec_x = elec_x 
        surf_elec_y = elec_y 
        if elec_z is None:
            surf_elec_z = np.zeros_like(elec_x)
            elec_z = np.zeros_like(elec_x)
        else:
            surf_elec_z = elec_z.copy()
            elec_z = np.array(elec_z) - np.array(elec_z)#normalise elec_z
        bur_elec_x = []
        bur_elec_y = []
        bur_elec_z = []
        x_interp = np.append(surf_elec_x,surf_x)
        y_interp = np.append(surf_elec_y,surf_y)
        z_interp = np.append(surf_elec_z,surf_z)
            
    #check if remeote electrodes present, and remove them 
    if len(rem_elec_idx)>0:
        elec_x = np.delete(elec_x,rem_elec_idx)
        elec_y = np.delete(elec_y,rem_elec_idx)
        elec_z = np.delete(elec_z,rem_elec_idx)  
    
    # check for some meshing parameters if there are aany 
    dist = np.unique(findDist(elec_x, elec_y, elec_z))
    if 'cl' in kwargs.keys(): 
        cl = kwargs['cl']
    else: 
        cl = dist[1]/4 
    
    if 'fmd' in kwargs.keys(): # compute depth of investigation if not given 
        fmd = kwargs['fmd']
    else:
        fmd = dist[-1]/3 # maximum possible dipole length / 3 (assuming a surface array)
    
    if 'cl_factor' in kwargs.keys():
        cl_factor = kwargs['cl_factor']
    else:
        cl_factor = 5 
        
    if 'cln_factor' in kwargs.keys():
        cln_factor = kwargs['cln_factor']
    else:
        cln_factor = 10000
        
    # check if mesh refinement present 
    if mesh_refinement is not None and not whole_space:        
        surf_rx = []
        surf_ry = []
        surf_rz = []
        surf_idx = []
        bur_rx = []
        bur_ry = []
        bur_rz = []
        bur_idx = []
        for i, key in enumerate(mesh_refinement['type']):
            if key == 'buried':
                bur_rx.append(mesh_refinement['x'][i])
                bur_ry.append(mesh_refinement['y'][i])
                bur_rz.append(mesh_refinement['z'][i])
                bur_idx.append(i)
            if key == 'surface' or key=='electrode':
                surf_rx.append(mesh_refinement['x'][i])
                surf_ry.append(mesh_refinement['y'][i])
                surf_rz.append(mesh_refinement['z'][i])
                surf_idx.append(i)
        #interpolate in order to normalise buried electrode elevations to 0
        x_interp = np.append(x_interp,surf_rx)#parameters to be interpolated with
        y_interp = np.append(y_interp,surf_ry)
        z_interp = np.append(z_interp,surf_rz)
        
        #if we have buried points normalise their elevation to as if they are on a flat surface
        if len(bur_rz)>0: 
            dump('found buried mesh refinement points')
            if interp_method is None:# still need to normalise electrode depths if we want a flat mesh, so use biliner interpolation instead
                bur_rz_topo = np.zeros_like(bur_idx)
            elif interp_method == 'bilinear': 
                bur_rz_topo = interp.interp2d(bur_rx, bur_ry, x_interp, y_interp, z_interp)
            elif interp_method == 'nearest':
                bur_rz_topo = interp.nearest(bur_rx, bur_ry, x_interp, y_interp, z_interp)
            elif interp_method == 'spline':
                bur_rz_topo = interp.interp2d(bur_rx, bur_ry, x_interp, y_interp, z_interp,method='spline')
            elif interp_method == 'triangulate':
                bur_rz_topo = interp.triangulate(bur_rx, bur_ry, x_interp, y_interp, z_interp)
                
        rz = np.array(mesh_refinement['z']) 
        rz[surf_idx] = 0
        if len(bur_rz)>0:
            rz[bur_idx] = rz[bur_idx] - bur_rz_topo # normalise to zero surface 
        
        internal_mesh_refinement = {'x':mesh_refinement['x'],
                                    'y':mesh_refinement['y'],
                                    'z':rz}
        if 'cl' in mesh_refinement.keys(): # pass individual characteristic lengths to box_3d 
            internal_mesh_refinement['cl'] = mesh_refinement['cl']
    
        kwargs['mesh_refinement'] = internal_mesh_refinement # actual mesh refinement passed to halfspace problem 
    elif mesh_refinement is not None and whole_space: 
        internal_mesh_refinement = {'x':mesh_refinement['x'],
                                    'y':mesh_refinement['y'],
                                    'z':mesh_refinement['z']}
        kwargs['mesh_refinement'] = internal_mesh_refinement # actual mesh refinement passed to wholespace problem
    elif whole_space:
        tree = cKDTree(np.c_[elec_x,elec_y,elec_z]) 
        cpx = np.zeros(len(bur_elec_x), dtype = float)
        cpy = np.zeros(len(bur_elec_x), dtype = float)
        cpz = np.zeros(len(bur_elec_x), dtype = float)
        # setup some code to get control points spiraling round the electrodes 
        r = cl*1.0
        choicex = [r, 0, -r, 0]
        choicey = [0, -r, 0, r]
        j = 0 
        for i in range(len(elec_x)):
            cpx[i] = elec_x[i] + choicex[j]
            cpy[i] = elec_y[i] + choicey[j]
            cpz[i] = elec_z[i] 
            j += 1 
            if j >= 4:
                j = 0 
        
        cpl = np.full_like(cpx, cl*cl_factor)
                
        # check that no points are duplicates of electrodes 
        idist,_ = tree.query(np.c_[cpx,cpy,cpz]) 
        tokeep = [True]*len(cpx)
        for i in range(len(cpx)):
            if idist[i] < 1e-16:
                tokeep[i] = False 
        cpx = cpx[tokeep]
        cpy = cpy[tokeep]
        cpz = cpz[tokeep]
        
        control = {'x':cpx,
                   'y':cpy,
                   'z':cpz,
                   'cl':cpl}
        kwargs['mesh_refinement'] = control 
    elif not whole_space: # make some of our own control points in this case
        elec = np.c_[elec_x, elec_y, elec_z] # all electrodes 
        tree = cKDTree(elec) 
        idist,_ = tree.query(elec,5) 
        
        cpx = np.array([])
        cpy = np.array([])
        cpz = np.array([])
        cpl = np.array([])
        
        ## control point parameters 
        cps = np.mean(idist[:,1]) # control point spacing 
        check_quadrant = True 
        # reduce cps for tightly spaced electrodes 
        tight_count = 0 
        total_count = 0 
        for i in range(idist.shape[0]):
            for j in range(1,idist.shape[1]):
                if abs(cps - idist[i,j]) < (cps*0.2):
                    tight_count += 1 
                total_count += 1 
        
        if (tight_count/total_count)>0.8:
            cps = cps/2
            # print('Tightly packed electrodes detected! Changing control point constraints')
            check_quadrant = False 
            
        if len(surf_elec_x)>0 and not coplaner: 
            # round 0 -> control points in x y coordinates surrounding the electrodes 
            _cpx, _cpy, _cpl = halfspaceControlPoints(surf_elec_x, surf_elec_y, 
                                                      cps, cl, check_quadrant, 
                                                      cl_factor, cln_factor)
            _cpz = np.zeros_like(_cpx)
            
            fmdx = _cpx[::3]
            fmdy = _cpy[::3]
            fmdz = _cpz[::3] - abs(fmd) # add some points at depth too 
            fmdl = _cpl[::3] 
        
            # append to control points 
            cpx = np.append(_cpx,fmdx)
            cpy = np.append(_cpy,fmdy)
            cpz = np.append(_cpz,fmdz)
            cpl = np.append(_cpl,fmdl)
        
        if len(bur_elec_x) > 0:
            # add some random control points near the electrodes 
            cbx = [0]*len(bur_elec_x)
            cby = [0]*len(bur_elec_x)
            cbz = [0]*len(bur_elec_x)
            cbl = [0]*len(bur_elec_x)
            # setup some code to get control points spiraling round the borehole electrodes 
            r = cl*1.0
            choicex = [r, 0, -r, 0]
            choicey = [0, -r, 0, r]
            j = 0 
            for i in range(len(bur_elec_idx)):
                cbx[i] = bur_elec_x[i] + choicex[j]
                cby[i] = bur_elec_y[i] + choicey[j]
                cbz[i] = bur_elec_z[i] 
                cbl[i] = cl*cl_factor 
                j += 1 
                if j >= 4:
                    j = 0 
                
            cpx = np.append(cpx,cbx)
            cpy = np.append(cpy,cby)
            cpz = np.append(cpz,cbz)
            cpl = np.append(cpl,cbl)
            # check that if electodes are close to the surface then add some 
            # control at the surface 
            min_dist_from_sur = abs(max(bur_elec_z))

            if min_dist_from_sur <= (cl*cl_factor):
                ux = np.unique(bur_elec_x)
                uy = np.unique(bur_elec_y)
                for i in range(len(ux)):
                    for j in range(len(uy)):
                        cpx = np.append(cpx,ux[i])
                        cpy = np.append(cpy,uy[j])
                        cpz = np.append(cpz,0)
                        cpl = np.append(cpl,cl)
                
        # one last check that no points are duplicates of electrodes 
        idist,_ = tree.query(np.c_[cpx,cpy,cpz]) 
        tokeep = [True]*len(cpx)
        for i in range(len(cpx)):
            if cpz[i] == 0 and len(surf_elec_x) == 0:
                continue 
            if idist[i] < (cps*0.8):
                tokeep[i] = False 
        cpx = cpx[tokeep]
        cpy = cpy[tokeep]
        cpz = cpz[tokeep]
        cpl = cpl[tokeep]
        
        control = {'x':cpx,
                   'y':cpy,
                   'z':cpz,
                   'cl':cpl}
        kwargs['mesh_refinement'] = control 
    
    
    # Use ball refinement by (default)
    if ball_refinement:
        kwargs['use_fields'] = True 
        
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
        node_pos = gw.wholespace3d(elec_x, elec_y, elec_z, file_path=file_name, **kwargs)
    elif lineis2d:
        node_pos = gw.halfspace3dline2d(elec_x, elec_y, elec_z, file_path=file_name, **kwargs)
    else:
        node_pos = gw.halfspace3d(elec_x, elec_y, elec_z, file_path=file_name, **kwargs)
            
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
    mesh.fmd = fmd 
    
    #mesh.write_dat(file_path='mesh.dat') # write mesh.dat - disabled as handled higher up in the R2 class 
    node_x = np.array(mesh.node[:,0])
    node_y = np.array(mesh.node[:,1])
    
    if keep_files is False: 
        os.remove(file_name+".geo");os.remove(file_name+".msh")
    
    if model_err:
        interp_method = None 
        dump('Created flat mesh for forward error modelling')
    else:
        dump('interpolating topography onto mesh using %s interpolation...'%interp_method)
    
    #using home grown functions to interpolate / extrapolate topography on mesh
    if interp_method == 'bilinear':# interpolate on a irregular grid, extrapolates the unknown coordinates
        nodez = interp.interp2d(node_x, node_y, x_interp, y_interp, z_interp)
    elif interp_method == 'nearest':
        nodez = interp.nearest(node_x, node_y, x_interp, y_interp, z_interp)
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
    mesh.elec_type = elec_type
    
    return mesh

#%% voxel mesh 
def voxelMesh(elec_x, elec_y, elec_z=None, elec_type = None, keep_files=True, 
              interp_method = 'triangulate', surface_refinement=None, 
              mesh_refinement=None, ball_refinement=True,
              path='exe', dump=print, whole_space=False, model_err=False,
              handle=None, show_output=True, **kwargs):
    """Generates a voxel mesh, not meant to be used for ERT processing with the R*
    codes but rather can be used for the purposes of mesh visualizuation post 
    processing, and is more readily exportable into formats required by
    other geological software (e.g. Geovisionary). 
    
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
        Not used. Kept as an argument for compatiblity with TetraMesh. 
    interp_method: string, default ='triangulate' optional
        Interpolation method to translate mesh nodes in the z direction. In 
        other words the method in which topography is appended to the mesh. 
        Here the topography is added to the mesh in post processing. 
        The options are documented below in the notes. 
    surface_refinement : np.array, optional 
        Numpy array of shape (3,n), should follow the format 
        np.array([x1,x2,x3,...],[y1,y2,y3,...],[z1,z2,z3,...]).
        Allows for extra refinement for the top surface of the mesh. 
        The points are not added to the mesh, but considered in post processing 
        in order to super impose topography on the mesh. 
    mesh_refinement : dict, pd.DataFrame, optional 
        Not used. Kept as an argument for compatiblity with TetraMesh. 
    ball_refinement: boolean, optional
        Not used. Kept as an argument for compatiblity with TetraMesh. 
    path : string, optional
        Not used. Kept as an argument for compatiblity with TetraMesh. 
    whole_space: boolean, optional
        Flag for if the problem should be treated as a whole space porblem, in 
        which case electrode type is ignored. 
    model_err: bool
        If True, a flat mesh will be returned for the sake of estimating 
        forward modelling errors. 
    dump : function, optional
        Function to which pass the output during mesh generation. `print()` is
        the default.
    **kwargs: dict
        Keyword arguments to be passed to functions in gmshWrap.py. Pass cl=x
        to force the voxel mesh characteristic length. 
    force_regular: bool, optional 
        If passed as True, then return voxel mesh will be made of elements which 
        are of the same size. 
         
    Returns
    -------
    mesh3d: class
    
    
    Notes 
    -----
    Possible arguments for interp_method: 
        'bilinear' : 4 known points are used to compute the equation of a plane 
                    in which the interpolated point lies. This method is reccomended 
                    if elevation data is organised in a regular grid. 
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
                    
    Format of mesh_refinement:
        The variable must have x,y,z and type columns
        'x': x coordinate array like 
        'y': y coordinate array like 
        'z': z coordinate array like 
        'type': list, object array of point types, use the tag 'surface' for 
            surface points, use 'buried' for buried points. 
        'cl': array like of characteristic lengths for each refinement point in 
            the mesh.
    """
    #formalities 
    if elec_z is None: 
        elec_z = np.zeros_like(elec_x)
        
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
            
    avail_methods = ['bilinear','nearest','spline','triangulate',None]
    if interp_method not in avail_methods:
        raise NameError("'%s' is an unrecognised interpretation method"%interp_method)

    if all(np.array(elec_z)==0):
        interp_method = None 
        
    if whole_space:
        warnings.warn('Voxel mesh is not yet fully tested for whole space problems')
        interp_method = None 
        
    if surface_refinement is not None:
        surf_x = surface_refinement[:,0]
        surf_y = surface_refinement[:,1]
        surf_z = surface_refinement[:,2]
    else:
        surf_x = []
        surf_y = []
        surf_z = []
    
    # setup electrodes 
    if elec_z is None:
        elec_z = np.zeros_like(elec_x)
        
    # check for repeated electrodes? 
    check4repeatNodes(elec_x,elec_y,elec_z,elec_type)
        
    elec = np.c_[elec_x,elec_y,elec_z]
    elec_cache = elec.copy() 
    
    rem_elec_idx = []
    if elec_type is not None and not whole_space:
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
                elec_cache[i,:] = -9999
    else: 
        surf_elec_x = elec_x 
        surf_elec_y = elec_y 
        surf_elec_z = elec_z 
        bur_elec_x = []
        bur_elec_y = []
        bur_elec_z = []
        
    #check if remeote electrodes present, and remove them 
    if len(rem_elec_idx)>0:
        elec_x = np.delete(elec_x,rem_elec_idx)
        elec_y = np.delete(elec_y,rem_elec_idx)
        elec_z = np.delete(elec_z,rem_elec_idx)  
        elec = np.c_[elec_x,elec_y,elec_z]

    #interpolate in order to normalise buried electrode elevations to 0
    x_interp = np.append(surf_elec_x,surf_x)#parameters to be interpolated with
    y_interp = np.append(surf_elec_y,surf_y)
    z_interp = np.append(surf_elec_z,surf_z)
    
    if len(x_interp) == 0:
        raise Exception('Cannot make a voxel mesh without any surface control points!') 
    
    if interp_method == 'triangulate':
        # need to check the number of interpolation points is stable for triangulation 
        if len(x_interp) == 4: 
            interp_method = 'bilinear'
        elif len(x_interp) < 4: 
            interp_method = 'nearest'
    
    if all(y_interp==y_interp[0]): # cant interpolate if on one plane 
        interp_method = 'nearest'
    
    if interp_method is None:
        elec_z_topo = np.zeros_like(elec_x) # still need to normalise electrode depths if we want a flat mesh, so use biliner interpolation instead
    elif interp_method == 'bilinear': 
        elec_z_topo = interp.interp2d(elec_x, elec_y, x_interp, y_interp, z_interp)
    elif interp_method == 'nearest':
        elec_z_topo = interp.nearest(elec_x, elec_y, x_interp, y_interp, z_interp)
    elif interp_method == 'spline':
        elec_z_topo = interp.interp2d(elec_x, elec_y, x_interp, y_interp, z_interp, method='spline')
    elif interp_method == 'triangulate':
        elec_z_topo = interp.triangulate(elec_x, elec_y, x_interp, y_interp, z_interp)
        
    ## mesh param 
    dp_len = np.max(findDist(elec_x, elec_y, elec_z))
    tree = cKDTree(elec) 
    idist,_ = tree.query(elec,2) 
    
    if 'cl' in kwargs.keys(): 
        cl = kwargs['cl']
    else: 
        cl = min(idist[:,1])/4 
    
    if 'fmd' in kwargs.keys() and kwargs['fmd'] is not None: # compute depth of investigation if not given 
        fmd = kwargs['fmd']
    else:
        fmd = dp_len/3 # maximum possible dipole length / 3
        if whole_space: 
            fmd = abs(max(elec_z) - min(elec_z))*1.1 
    
    if 'cl_factor' in kwargs.keys() and kwargs['cl_factor'] is not None:
        cl_factor = kwargs['cl_factor']
    else:
        cl_factor = 5 
        
    if 'force_regular' in kwargs.keys():
        force_regular = kwargs['force_regular']
    else:
        force_regular = False 

    elec_z = elec_z - elec_z_topo
    elec_z = np.round(elec_z,6) # nuke near zero values 
    
    ## node generation algorithm 
    dump('Running voxel mesh node generation...')
    # generate x - y plane first, will add topo later 
    ex = np.unique(elec_x).tolist()
    ey = np.unique(elec_y).tolist()
    ez = np.unique(elec_z).tolist()
    ex = [ex[0] - cl*cl_factor] + ex 
    ex = ex + [ex[-1] + cl*cl_factor]
    ey = [ey[0] - cl*cl_factor] + ey 
    ey = ey + [ey[-1] + cl*cl_factor]
    ez = [max(ez) - fmd] + ez 
    
    ## insert x values 
    xx = []
    x = min(ex)
    for i in range(1,len(ex)):
        while ex[i] > x:
            x += cl 
            xx.append(x)
        if not force_regular: 
            xx.append(ex[i])
        
    dump('Generated %i unique nodes in X direction'%len(xx))
    
    ## insert y values 
    yy = []
    y = min(ey)
    for i in range(1,len(ey)):
        while ey[i] > y:
            y += cl 
            yy.append(y)
        if not force_regular: 
            yy.append(ey[i])
    
    dump('Generated %i unique nodes in Y direction'%len(yy))
        
    ## insert z values 
    zz = [0]
    z = min(ez)
    for i in range(1,len(ez)):
        while ez[i] > z and (z + cl) <=0:
            z += cl 
            zz.append(z)
        if not force_regular: 
            zz.append(ez[i])
        
    dump('Generated %i unique nodes in Z direction'%len(zz)) 
    
    ## force unique sorting 
    xx = np.unique(xx)
    yy = np.unique(yy)
    zz = np.unique(zz)

    # meshx, meshy, meshz = np.meshgrid(xx,yy,zz) # old line --> works 
    
    # changed node creation to be friendly to vts format 
    meshz, meshy, meshx = np.meshgrid(zz,yy,xx,indexing='ij')
    # need to check this does not mess up the connection matrix creation... ?
    
    # compress node arrays into columns 
    node = np.c_[meshx.flatten(),meshy.flatten(),meshz.flatten()]
    numnp = node.shape[0]

    ## add topography back to nodes 
    if model_err: # dont add any topography to the mesh in the case of modelling errors 
        interp_method = None 
        
    #using home grown functions to interpolate / extrapolate topography on mesh
    if interp_method == 'bilinear':# interpolate on a irregular grid, extrapolates the unknown coordinates
        mesh_z_topo = interp.interp2d(node[:,0], node[:,1], x_interp, y_interp, z_interp)
    elif interp_method == 'nearest':
        mesh_z_topo = interp.nearest(node[:,0], node[:,1], x_interp, y_interp, z_interp)
    elif interp_method == 'spline':
        mesh_z_topo = interp.interp2d(node[:,0], node[:,1], x_interp, y_interp, z_interp,method='spline')
    elif interp_method == 'triangulate':
        mesh_z_topo = interp.triangulate(node[:,0], node[:,1], x_interp, y_interp, z_interp)
    elif interp_method == None:
        mesh_z_topo = np.zeros_like(node[:,0],dtype=float)
    node[:,2] += mesh_z_topo 
    
    ## generate elements 
    # map the 8 nodes surrounding every element 
    nid = np.arange(numnp)
    numel = (len(xx)-1)*(len(yy)-1)*(len(zz)-1)
    nmesh = nid.reshape(meshx.shape) # node mesh 
    connec = np.zeros((numel,8),dtype=int)-1
    
    ii = 0 
    for i in range(nmesh.shape[0]-1):
        for j in range(nmesh.shape[1]-1):
            for k in range(nmesh.shape[2]-1):
                connec[ii,0] = nmesh[i,j,k]
                connec[ii,1] = nmesh[i+1,j,k]
                connec[ii,2] = nmesh[i,j+1,k]
                connec[ii,3] = nmesh[i+1,j+1,k]
                connec[ii,4] = nmesh[i,j,k+1]
                connec[ii,5] = nmesh[i+1,j,k+1]
                connec[ii,6] = nmesh[i,j+1,k+1]
                connec[ii,7] = nmesh[i+1,j+1,k+1]
                ii += 1 
                
    dump('Mesh consists of %i nodes and %i elements'%(numnp,numel)) 
    
    ## map electrodes to nodes 
    tree = cKDTree(node)
    _, enodes = tree.query(elec_cache)
    
    ## make mesh class    
    mesh = Mesh(node[:,0], node[:,1], node[:,2], 
                connec,
                cell_type=[11],
                original_file_path='N/A',
                order_nodes=False) 
    
    # add electrode nodes 
    iremote = np.array([False]*elec_cache.shape[0],dtype=bool) 
    for i in rem_elec_idx:
        iremote[i] = True 
    
    mesh.setElecNode(enodes, iremote)
    mesh.elec_type = elec_type
    
    #set the number of voxels in the xyz directions 
    mesh.nvoxel['x'] = len(xx)-1 
    mesh.nvoxel['y'] = len(yy)-1 
    mesh.nvoxel['z'] = len(zz)-1 
    
    return mesh 


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
    gw.prism([elec_x,elec_y,elec_z], file_path=file_name, **kwargs)
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
    

#%% cylinder mesh
def cylinderMesh(elec_x, elec_y, elec_z,
                 zlim=None, radius=None, file_path='cylinder_mesh.geo',
                 cl=-1, cl_factor=2, finer=4,
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
        circle shape. No longer used. 
    handle : variable, optional
        Will be assigned the output of 'Popen' in case the process needs to be
        killed in the UI for instance.
    """
    file_path = file_path if file_path[-4:] == '.geo' else file_path + '.geo'
    elec = [elec_x,elec_y,elec_z]
    gw.cylinder(elec, zlim=zlim, radius=radius, file_path=file_path,
                cl=cl, cl_factor=cl_factor, elemz=finer)    
    elec = np.array(elec).T 
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
# mesh = cylinderMesh(elec, file_path='invdir/mesh_cylinder',
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
    gw.tank(elec=elec, origin=origin, dimension=dimension,
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
    elif ext == '.vtu':
        mesh = vtuImport(file_path, order_nodes=order_nodes)
    elif ext == '.vts':
        mesh = vtsImport(file_path, order_nodes=False)
    elif ext == '.msh':
        mesh_dict = gw.mshParse(file_path, debug=False)

        mesh = Mesh(node_x = mesh_dict['node_x'],
                    node_y = mesh_dict['node_y'],
                    node_z = mesh_dict['node_z'],
                    node_data = np.array(mesh_dict['node_data']).T,
                    cell_type = mesh_dict['cell_type'],
                    order_nodes = order_nodes)
        mesh.addAttribute(mesh_dict['parameters'], 'zones')
        mesh.addAttribute(mesh_dict['parameters'], 'region')
        
    elif ext == '.dat':
        mesh = dat_import(file_path, order_nodes=order_nodes)   
    elif ext == '.node':
        mesh = tetgen_import(file_path, order_nodes=order_nodes)
    else:
        avail_ext = ['.vtk','.msh','.dat','.node']
        raise ImportError("Unrecognised file extension, available extensions are "+str(avail_ext))
    
    if node_pos is not None:
        mesh.setElecNode(np.array(node_pos,dtype=mesh.dint)) # add electrode nodes to mesh provided by the user
    
    # check if poorly formed (will raise an error if so)
    flag = ['mesh node']*mesh.numnp 
    check4repeatNodes(mesh.node[:,0], mesh.node[:,1], mesh.node[:,2], flag)
    
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

#%% write descrete points to a vtk file 
def points2vtk (x,y,z,fname="points.vtk",title='points',data=None, file_name=None):
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
    fname : string, optional
        Path to saved file, defualts to 'points.vtk' in current working directory.
    title : string, optional
        Title of vtk file.
    data: dict, pd.DataFrame, optional 
        Point data. 
    file_name: str, optional
        Same input as fname but kept for backwards compatibility 
            
    Returns
    -------
    ~.vtk : file
    """
    #error check
    if len(x) != len(y) or len(x) != len(z):
        raise ValueError('mis-match between vector lengths')
        
    if file_name is not None:
        fname = file_name 
    
    fh=open(fname,'w');#open file handle
    #add header information
    fh.write('# vtk DataFile Version 3.0\n')
    fh.write(title+'\n')
    fh.write('ASCII\n')
    fh.write('DATASET POLYDATA\n')
    #add data
    fh.write('POINTS      %i double\n'%len(x))
    [fh.write('{:<10} {:<10} {:<10}\n'.format(x[i],y[i],z[i])) for i in range(len(x))]
    
    if data is not None and len(data.keys())>0:
        fh.write("POINT_DATA %i\n"%len(x))
        for i,key in enumerate(data.keys()):
            fh.write("SCALARS %s double 1\n"%key.replace(' ','_'))
            fh.write("LOOKUP_TABLE default\n")
            X = np.array(data[key])
            X[np.isnan(X)]=-9999
            [fh.write("%16.8f "%X[j]) for j in range(len(x))]
            fh.write("\n")
    
    fh.close() 
    
def points2vtp (x,y,z,fname="points.vtp",title='points',data=None,file_name=None, 
                connect=True):
    """
    Function makes a .vtp file for some xyz coordinates. optional argument
    renames the name of the file (needs file path also) (default is "points.vtp"). 
            
    Parameters
    ----------
    x : list, tuple, np array
        X coordinates of points.
    y : list, tuple, np array
        Y coordinates of points.
    z : list, tuple, np array
        Z coordinates of points.
    fname : string, optional
        Path to saved file, defualts to 'points.vtk' in current working directory.
    title : string, optional
        Title of vtk file.
    data: dict, pd.DataFrame, optional 
        Point data. 
    file_name: str, optional
        Same input as fname but kept for backwards compatibility 
    connect: bool, optional 
        Connect the electrodes in a line. Can be helpful for paraview visualisation 
        options. 
            
    Returns
    -------
    ~.vtk : file
    """
    #error check
    if len(x) != len(y) or len(x) != len(z):
        raise ValueError('mis-match between vector lengths')
    
    if file_name is not None:
        fname = file_name 
    
    fh = open(fname,'w')
    
    def writeXMLline(text,bracket=True,tab=0):
        for i in range(tab):
            fh.write('\t')
        if bracket: 
            fh.write('<')
        fh.write('{:}'.format(text))
        if bracket:
            fh.write('>')
        fh.write('\n')
        
    def writeXMLarray(X):
        # tab = 5
        j = 0 
        text = '' 
        for i in range(len(X)): 
            if isinstance(X[i],int):
                text += '{:d} '.format(X[i])
            else: 
                text += '{:f} '.format(X[i])
            j += 1 
            if j == 5: 
                writeXMLline(text,False,5)
                j = 0 
                text = ''  
        if j < 5: # write out last line 
            writeXMLline(text,False,5)
            
    def writeEmpty(tag): 
        writeXMLline(tag,tab=3)
        writeXMLline(empty_header_template.format('connectivity'),tab=4)
        writeXMLline('/DataArray',tab=4)
        writeXMLline(empty_header_template.format('offsets'),tab=4)
        writeXMLline('/DataArray',tab=4)
        writeXMLline('/'+tag,tab=3)
    
    data_header_template = 'DataArray type="Float64" Name="{:s}" format="ascii" RangeMin="{:f}" RangeMax="{:f}"'
    data_header_templatei = 'DataArray type="Int64" Name="{:s}" format="ascii" RangeMin="{:d}" RangeMax="{:d}"'
    empty_header_template = 'DataArray type="Int64" Name="{:s}" format="ascii" RangeMin="1e+299" RangeMax="-1e+299"'
    point_header_template = 'DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii" RangeMin="{:f}" RangeMax="{:f}"'
   
    writeXMLline('VTKFile type="PolyData" version="1.0" byte_order="LittleEndian" header_type="UInt64"')
    writeXMLline('PolyData',tab=1)
    nlines = 0 
    if connect:
        nlines = 1 
    writeXMLline('Piece NumberOfPoints="%i" NumberOfVerts="0" NumberOfLines="%i" NumberOfStrips="0" NumberOfPolys="0"'%(len(x),nlines), tab=2)
    
    # write out point data
    writeXMLline('PointData', tab=3)
    if data is not None: 
        for column in data.keys():
            X = data[column] 
            header = data_header_template.format(column, np.min(X), np.max(X))
            writeXMLline(header, tab=4)
            X[np.isnan(X)] = -9999
            writeXMLarray(X)
            writeXMLline('/DataArray',tab=4)
    writeXMLline('/PointData', tab=3)
    
    # write cell data 
    writeXMLline('CellData',tab=3)
    writeXMLline('/CellData',tab=3)
    
    # write out node/point coordinates 
    writeXMLline('Points',tab=3)
    points = np.c_[x,y,z]
    header = point_header_template.format(np.min(points.flatten()), np.max(points.flatten()))
    writeXMLline(header, tab=4)
    for i in range(len(x)):
        line = '{:f}\t{:f}\t{:f}'.format(points[i,0], points[i,1], points[i,2])
        writeXMLline(line, False, 5)
    writeXMLline('/DataArray',tab=4)
    writeXMLline('/Points',tab=3)
    
    # write out vertices (not applicable here)
    writeEmpty('Verts')
    
    # connect up the electrodes if requested. 
    if connect: 
        writeXMLline('Lines', tab=3)
        X = [i for i in range(len(x))]
        header = data_header_templatei.format("connectivity", np.min(X), np.max(X))
        writeXMLline(header, tab=4)
        writeXMLarray(X)
        writeXMLline('/DataArray',tab=4)
        writeXMLline(data_header_templatei.format("offsets", len(x), len(x)),tab=4)
        writeXMLline('%i'%len(x),False,5)
        writeXMLline('/DataArray',tab=4)
        writeXMLline('/Lines',tab=3)
    else: 
        writeEmpty('Lines')
    
    # write out other data arrays which are empty 
    writeEmpty('Strips')
    writeEmpty('Polys')

    
    # closing lines 
    writeXMLline('/Piece',tab=2)
    writeXMLline('/PolyData',tab=1)
    writeXMLline('/VTKFile')
    
    fh.close()
    
    
    
