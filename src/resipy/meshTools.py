# MESH TOOLS 
"""
Created on Wed May 30 10:19:09 2018, python 3.6.5
@author: jamyd91
Module handles mesh generation, display, discretisation and post processing. 
The convention for x y z coordinates is that the z coordinate is the elevation.

Dependencies: 
    numpy (conda lib)
    matplotlib (conda lib)
    gmshWrap(pyR2 resipy module)
    python3 standard libaries
"""
#import standard python packages
import os, platform, warnings, multiprocessing, re, pathlib
from subprocess import PIPE, Popen, call
import time
#import matplotlib and numpy packages 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import ListedColormap
import matplotlib.tri as tri
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
#import R2gui API packages - remove the "api." in below code if running the script and using the test blocks at the bottom 
import resipy.gmshWrap as gw
from resipy.isinpolygon import isinpolygon, isinvolume, in_box
import resipy.interpolation as interp
from resipy.sliceMesh import sliceMesh # mesh slicing function

        
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
    cax = None 
    zone = None
    attr_cache={}
    mesh_title = "2D_R2_mesh"
    no_attributes = 0
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
                 elm_area,#area of each element (actually this can be left blank)
                 cell_type,#according to vtk format
                 cell_attributes,#the values of the attributes given to each cell 
                 atribute_title,#what is the attribute? we may use conductivity instead of resistivity for example
                 original_file_path='N/A',
                 regions=None) :

        #assign varaibles to the mesh object 
        self.num_nodes=num_nodes
        self.num_elms=num_elms
        self.node_x = node_x
        self.node_y = node_y
        self.node_z = node_z
        self.node_id=node_id
        self.elm_id=elm_id
        self.con_matrix = node_data #connection matrix
        self.elm_centre=elm_centre
        self.elm_area=elm_area
        self.cell_type=cell_type
        self.cell_attributes=cell_attributes 
        self.atribute_title=atribute_title
        self.original_file_path=original_file_path
        self.regions = regions
        self.surface = None # surface points for cropping the mesh when contouring
        #decide if mesh is 3D or not 
        if max(node_y) - min(node_y) == 0: # mesh is probably 2D 
            self.ndims=2
        else:
            self.ndims=3
            self.mesh_title = '3D_R3t_mesh' 
    
    @classmethod # creates a mesh object from a mesh dictionary
    def mesh_dict2class(cls, mesh_info):
        """ Converts a mesh dictionary produced by the gmsh2r2mesh and
        vtk_import functions into a mesh object, its an alternative way to
        make a mesh object. 
        ***Intended for development use***
            
        Parameters
        ----------
        mesh_info: dictionary 
            mesh parameters stored in a dictionary rather than a mesh, useful for debugging parsers
            
        Returns
        -------
        Mesh: class 
        """
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
                     mesh_info['original_file_path'])
        try:
            obj.add_attr_dict(mesh_info['cell_attributes'])
        except KeyError as e:
            #print('error in add_attr_dict', e)
            pass
        try:
            obj.regions = mesh_info['element_ranges']
        except KeyError:
            pass
                
        return (obj)
    

    def add_e_nodes(self,e_nodes):
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
        
    def file_path(self):#returns the file path from where the mesh was imported
        return(format(self.original_file_path))
       
    def Type2VertsNo(self):#converts vtk cell types into number of vertices each element has 
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
        out += "Number of cell vertices: %i\n"%self.Type2VertsNo()
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
#            print(len(values),self.num_elms)
            raise ValueError("The length of the new attributes array does not match the number of elements in the mesh")
        self.no_attributes += 1
        try: 
            self.attr_cache[key]=values #allows us to add an attributes to each element.
        except AttributeError:
            self.attr_cache = {}
            self.attr_cache[key]=values #add attribute 
    
    def add_attr_dict(self,attr_dict):
        """Mesh attributes are stored inside a dictionary, mesh.attr_cache.
        
        Parameters
        ------------
        attr_dict: dict 
            Each key in the dictionary should reference an array like of values. 
        """
        self.attr_cache=attr_dict
        self.no_attributes = len(attr_dict)
        
    def show_avail_attr(self,flag=True):
        """Show available attributes in mesh.attr_cache. 
        """
        out = '\n______cell attributes_____\n'
        try: 
            for i,key in enumerate(self.attr_cache):
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
    

    def show(self,color_map = 'Spectral',#displays the mesh using matplotlib
             color_bar = True,
             xlim = "default",
             zlim = "default",
             ax = None,
             electrodes = True,
             sens = False,
             edge_color = 'k',
             contour=False,
             vmin=None,
             vmax=None,
             attr=None,
             clabel=None,
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
        clabel: string, optional
            Label of the colorbar. Default is the value of `attr` argument.

        Returns
        ----------
        figure : matplotlib figure 
            Figure handle for the plotted mesh object.
        
        Notes
        ----------
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
            
        if self.ndims == 3:
            self.show_3D(color_map = color_map,#displays the mesh using matplotlib
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
            
        
        #decide which attribute to plot, we may decide to have other attritbutes! 
        if attr is None: 
            #plots default attribute
            X=np.array(self.cell_attributes) # maps resistivity values on the color map
            color_bar_title = self.atribute_title
        else:
            try:
                X = np.array(self.attr_cache[attr])
                color_bar_title = attr
            except (KeyError, AttributeError):
                raise KeyError("Cannot find attr_cache attribute in mesh object or 'attr' does not exist.")
        
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
            if xlim=="default":
                xlim=[min(self.elec_x),max(self.elec_x)]
            if zlim=="default":
                doiEstimate = 2/3*np.abs(self.elec_x[0]-self.elec_x[-1]) # TODO depends on longest dipole
                #print(doiEstimate)
                zlim=[min(self.elec_z)-doiEstimate,max(self.elec_z)]
        except AttributeError:
            if xlim=="default":
                xlim=[min(self.node_x),max(self.node_x)]
            if zlim=="default":
                zlim=[min(self.node_z),max(self.node_z)]
                
        ##plot mesh! ##
        a = time.time() #start timer on how long it takes to plot the mesh
        #compile mesh coordinates into polygon coordinates  
        nodes = np.c_[self.node_x, self.node_z]
        connection = np.array(self.con_matrix).T # connection matrix 
        #compile polygons patches into a "patch collection"
        ###X=np.array(self.cell_attributes) # maps resistivity values on the color map### <-- disabled 
        coordinates = nodes[connection]
        if vmin is None:
            vmin = np.min(X)
        if vmax is None:
            vmax = np.max(X)
#        if attr is not None:
#            if 'difference' in attr or 'Difference' in attr:
#                vext = np.max([np.abs(vmin), np.abs(vmax)])
#                vmin = -vext
#                vmax = vext
#                color_map = 'bwr'
        
        if edge_color == None or edge_color=='none' or edge_color=='None':
            edge_color='face'#set the edge colours to the colours of the polygon patches

        if contour is False:
            if attr is None: # so the default material
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
            xc = np.array(self.elm_centre[0])
            yc = np.array(self.elm_centre[2])
            zc = np.array(X)
            x = np.array(self.node_x)
            y = np.array(self.node_z)
            
            # set scale arrangement
            if vmin is None:
                vmin = np.nanmin(zc)
            if vmax is None:
                vmax = np.nanmax(zc)
            if vmax > vmin:
                levels = np.linspace(vmin, vmax, 14) # to have 2 contours between two cbar values!
            else:
                levels = None
            
            if self.cell_type[0] == 9: # quadrilateral mesh (exact topo)   
                # interpolate the cell-centered value to the node to be able
                # to use the triangular mesh already in the grid
                z = interp.nearest(x, y, xc, yc, zc)
                
                def rebuildRegularGrid(x2, y2, z2):
                    x2unique = np.unique(x2)
                    xs = []
                    for xuni in x2unique: # for loop otherwise dataframe groupby
                        xs.append(np.where(x2 == xuni)[0])
                    minLength = np.min([len(a) for a in xs])
                    xs2 = []
                    for a in xs:
                        isort = np.argsort(y2[a])
                        xs2.append(a[isort][-minLength:])
                    xs2 = np.vstack(xs2)
                    X = x2[xs2]
                    Y = y2[xs2]
                    Z = z2[xs2]
                    return X, Y, Z
            
                Xi, Yi, Zi = rebuildRegularGrid(x, y, z)
                self.cax = ax.contourf(Xi, Yi, Zi, levels=levels, cmap=color_map)
            
            
            elif self.cell_type[0] == 5: # triangular mesh (exact topo)
                # interpolate the cell-centered value to the node to be able
                # to use the triangular mesh already in the grid
                z = interp.nearest(x, y, xc, yc, zc)
                triang = tri.Triangulation(x, y, connection)
                
                self.cax = ax.tricontourf(triang, z, levels=levels, extend='both')
            
            else: # fallback mode with tricontourf and cropSurface() (topo based on centroids) 
                triang = tri.Triangulation(xc, yc) # build grid based on centroids
                z = zc

                # make sure none of the triangle centroids are above the
                # line of electrodes
                def cropSurface(triang, xsurf, ysurf):
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
                    return i2keep
                
                try:
                    triang.set_mask(~cropSurface(triang, self.surface[:,0], self.surface[:,1]))
                except Exception as e:
                    print('Error in Mesh.show for contouring: ', e)
                
                self.cax = ax.tricontourf(triang, z, levels=levels, extend='both')
            
            
        ax.autoscale()
        #were dealing with patches and matplotlib isnt smart enough to know what the right limits are, hence set axis limits 
        ax.set_ylim(zlim)
        ax.set_xlim(xlim)
        ax.set_xlabel('Distance [m]')
        ax.set_ylabel('Elevation [m]')
        
        if color_bar:#add the color bar 
            self.cbar = plt.colorbar(self.cax, ax=ax, format='%.1f')
            self.cbar.set_label(color_bar_title) #set colorbar title

        ax.set_aspect('equal')#set aspect ratio equal (stops a funny looking mesh)

        #biuld alpha channel if we have sensitivities 
        if sens:
            try:
                weights = np.array(self.sensitivities) #values assigned to alpha channels 
                alphas = np.linspace(1, 0, self.num_elms)#array of alpha values 
                raw_alpha = np.ones((self.num_elms,4),dtype=float) #raw alpha values 
                raw_alpha[..., -1] = alphas
                alpha_map = ListedColormap(raw_alpha) # make a alpha color map which can be called by matplotlib
                #make alpha collection
                alpha_coll = PolyCollection(coordinates, array=weights, cmap=alpha_map, edgecolors='none', linewidths=0)#'face')
                #*** the above line can cuase issues "attribute error" no np.array has not attribute get_transform, 
                #*** i still cant figure out why this is becuase its the same code used to plot the resistivities 
                ax.add_collection(alpha_coll)
                
            except AttributeError:
                print("no sensitivities in mesh object to plot")
        
        if electrodes: #try add electrodes to figure if we have them 
            try: 
                ax.plot(self.elec_x,self.elec_z,'ko')
            except AttributeError:
                print("no electrodes in mesh object to plot")

        # adding interactive display when mouse-over
        centroids = np.array([self.elm_centre[0], self.elm_centre[2]]).T
        def format_coord(x, y):
            dist = np.sqrt(np.sum((centroids - np.array([x, y]))**2, axis=1))
            imin = np.argmin(dist)
            return ('x={:.2f} m, elevation={:.2f} m, value={:.3f}'.format(x,y,X[imin]))
        ax.format_coord = format_coord

        print('Mesh plotted in %6.5f seconds'%(time.time()-a))
        
        if iplot == True:
            return fig
    
    def draw(self, 
             attr=None,
             edge_color = 'k',
             color_map = None,
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
            
        if attr is None: 
            #plots default attribute
            X=np.array(self.cell_attributes) # maps resistivity values on the color map
            color_bar_title = self.atribute_title
        else:
            try:
                X = np.array(self.attr_cache[attr])
                color_bar_title = attr
            except (KeyError, AttributeError):
                raise KeyError("Cannot find attr_cache attribute in mesh object or 'attr' does not exist.")
                
        a = time.time() #start timer on how long it takes to plot the mesh
        
        if edge_color == None or edge_color=='none' or edge_color=='None':
            edge_color='face'#set the edge colours to the colours of the polygon patches
            
        if vmin is None:
            vmin = np.min(X)
        if vmax is None:
            vmax = np.max(X)
        
        if color_map != None :
            if attr is None:
                cm = plt.get_cmap(color_map, len(np.unique(X)))
                vmin = vmin + 0.5
                vmax = vmax - 0.5
            else:
                cm = color_map
            self.cax.set_cmap(cm) # change the color map if the user wants to 
        else:
            if attr is None:
                cm = plt.get_cmap('Spectral', len(np.unique(X)))
                vmin = vmin + 0.5
                vmax = vmax - 0.5
                self.cax.set_cmap(cm)
        
        
        #following block of code redraws figure 
        self.cax.set_array(X) # set the array of the polygon collection to the new attribute 
        self.cax.set_clim(vmin=vmin, vmax=vmax) # reset the maximum limits of the color map 
        self.ax.add_collection(self.cax)#blit polygons to axis
        self.cbar.set_label(color_bar_title) # change the color bar title 
        self.fig.canvas.draw() # redraw figure canvas (does not make a new figure it is faster fig.show())

        if color_bar:#add the color bar 
           print("you should have decided you wanted a color bar when using the mesh.show function")
            
        print('Mesh plotted in %6.5f seconds'%(time.time()-a))    
    
    
    def show_3D(self,color_map = 'Spectral',#displays the mesh using matplotlib
             color_bar = True,
             xlim = "default",
             ylim = "default",
             zlim = "default", 
             ax = None,
             electrodes = True,
             sens = False,
             edge_color = 'k',
             alpha = 1,
             vmax=None,
             vmin=None,
             attr=None):
        """
        Shows a 3D tetrahedral mesh. 
        
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

        Returns
        ----------
        figure : matplotlib figure 
            Figure handle for the plotted mesh object.
        
        Notes
        ----------
        Show a mesh object using matplotlib. The color map variable should be 
        a string refering to the color map you want (default is "jet").
        As we're using the matplotlib package here any color map avialable within 
        matplotlib package can be used to display the mesh here also. See: 
        https://matplotlib.org/2.0.2/examples/color/colormaps_reference.html
        
        Plotting sensitivies using sens=True is not reccomended. The matplotlib renderer has trouble with it. 
        """
        if not isinstance(color_map,str):#check the color map variable is a string
            raise NameError('color_map variable is not a string')
        
        if self.ndims==2:
            warnings.warn("Its reccomended to use mesh.show() for 2D meshes, results of 3D show will be unstable")
            
        #decide which attribute to plot, we may decide to have other attritbutes! 
        if attr is None: 
            #plots default attribute
            X=np.array(self.cell_attributes) # maps resistivity values on the color map
            color_bar_title = self.atribute_title
        else:
            try:
                X = np.array(self.attr_cache[attr])
                color_bar_title = attr
            except (KeyError, AttributeError):
                raise KeyError("Cannot find attr_cache attribute in mesh object or 'attr' does not exist.")
        
        t0 = time.time() # benchmark function
        
        #make 3D figure 
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        
        if xlim=="default":
            xlim=[min(self.elec_x), max(self.elec_x)]
        if ylim=="default":
            ylim=[min(self.elec_y), max(self.elec_y)]
        if zlim=="default":
            zlim=[min(self.node_z),max(self.node_z)]
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
        
        # search through each element to see if it is on the edge of the mesh, 
        # this step is important as it is very expensive to plot anything in 3D using matplotlib 
        # triangles on the edge of the mesh will be used only once
        tri_combo = np.zeros((self.num_elms,4),dtype='float64')
        elm_x = self.elm_centre[0]
        elm_y = self.elm_centre[1]
        elm_z = self.elm_centre[2]
        in_elem = in_box(elm_x,elm_y,elm_z,xlim[1],xlim[0],ylim[1],ylim[0],zlim[1],zlim[0])#find elements veiwable in axis
        X = X[in_elem]#reassign X to elements inside the box limits 
        temp_con_mat = np.array(self.con_matrix,dtype='int64')#temporary connection matrix which is just the elements inside the box
        con_mat=temp_con_mat[:,in_elem] # truncate elements
        inside_numel = len(con_mat[0])#number of inside elements 
        
        for i in range(inside_numel):
            idx1 = con_mat[0][i]#extract indexes 
            idx2 = con_mat[1][i]
            idx3 = con_mat[2][i]
            idx4 = con_mat[3][i]
            
            face1 = idx1*idx2*idx3 # hopefully each of these make a unique value 
            face2 = idx1*idx2*idx4
            face3 = idx2*idx3*idx4
            face4 = idx1*idx4*idx3

            tri_combo[i,0] = face1#face 1 
            tri_combo[i,1] = face2#face 2 
            tri_combo[i,2] = face3#face 3 
            tri_combo[i,3] = face4#face 4 
            
        #shape = tri_combo.shape
        tri_combo = tri_combo.flatten()
        temp,index,counts = np.unique(tri_combo,return_index=True,return_counts=True) # find the unique values 
        single_vals_idx = counts==1
        edge_element_idx = index[single_vals_idx]/4
        face_element_idx = np.floor(edge_element_idx)
        face_probe = edge_element_idx - np.floor(edge_element_idx)
        
        truncated_numel = len(face_element_idx)
        face_list = [0] * truncated_numel
        assign = [0] * truncated_numel # the number assigned to each face
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
            
            face1 = (vert1,vert2,vert3)
            face2 = (vert1,vert2,vert4)
            face3 = (vert2,vert3,vert4)
            face4 = (vert1,vert4,vert3)
            
            if face_probe[i] == 0: #if single_val_idx. == 0 > face1
                face_list[i] = face1#face 1 
            elif face_probe[i] == 0.25:#if single_val_idx. == 0.25 > face2
                face_list[i] = face2#face 2
            elif face_probe[i] == 0.5:#if single_val_idx. == 0.5 > face3
                face_list[i] = face3#face 3 
            elif face_probe[i] == 0.75:#if single_val_idx. == 0.75 > face4
                face_list[i] = face4#face 4
            
            assign[i] = X[ref]#get attribute value into assigned array
          
        polly = Poly3DCollection(face_list,linewidth=0.5) # make 3D polygon collection
        polly.set_alpha(alpha)#add some transparancy to the elements
        try:
            polly.set_array(np.array(assign))
        except MemoryError:#catch this error and print something more helpful than matplotlibs output
            raise MemoryError("Memory access voilation encountered when trying to plot mesh, \n please consider truncating the mesh or display the mesh using paraview.")
        polly.set_edgecolor(edge_color)
        polly.set_cmap(color_map) # set color map 
        polly.set_clim(vmin=vmin, vmax=vmax) # reset the maximum limits of the color map 
        ax.add_collection3d(polly, zs='z')#blit polygons to axis 
        self.cax = polly
        
        if color_bar:#add the color bar 
            self.cbar = plt.colorbar(self.cax, ax=ax, format='%.1f')
            self.cbar.set_label(color_bar_title) #set colorbar title
            
#        ax.set_aspect('equal')#set aspect ratio equal (stops a funny looking mesh)
        
        if sens: #add sensitivity to plot if available
            try:
                weights = np.array(self.sensitivities) #values assigned to alpha channels 
                alphas = np.linspace(1, 0, self.num_elms)#array of alpha values 
                raw_alpha = np.ones((self.num_elms,4),dtype=float) #raw alpha values 
                raw_alpha[..., -1] = alphas
                alpha_map = ListedColormap(raw_alpha) # make a alpha color map which can be called by matplotlib
                #make alpha collection
                alpha_coll = Poly3DCollection(face_list, array=weights, cmap=alpha_map, edgecolors='none', linewidths=0)#'face')
                ax.add_collection(alpha_coll)
            except AttributeError:
                print("no sensitivities in mesh object to plot")
            
        if electrodes: #try add electrodes to figure if we have them 
            try: 
                ax.scatter(self.elec_x,self.elec_y,zs=np.array(self.elec_z),
                           s=20, c='k', marker='o')#note you have to give the points a size otherwise you
                #get an NoneType Attribute error. 
                #the matplotlib renderer really doesn't cope well with the addition of the electrodes, 
                #the points are usually masked by the elements...
                #I reccomend putting an alpha setting on the mesh to view electrodes + mesh together. 
            except AttributeError as e:
                print("could not plot 3d electrodes, error = "+str(e))
            
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
        
        
    def assign_zone(self,poly_data):
        """ Assign material/region assocations with certain elements in the mesh 
        say if you have an area you'd like to forward model. 
        ***2D ONLY***
            
        Parameters
        ----------
        poly_data : dictionary 
            Dictionary with the vertices (x,y) of each point in the polygon.
            
        Returns
        ---------
        material_no : numpy.array
            Element associations starting at 1. So 1 for the first region 
            defined in the region_data variable, 2 for the second region 
            defined and so on. If the element can't be assigned to a region
            then it'll be left at 0. 
        """   
        no_elms=self.num_elms#number of elements 
        elm_xyz=self.elm_centre#centriods of mesh elements 
        material_no=np.zeros(no_elms,dtype=int)#attribute number
        
        if not isinstance(poly_data,dict):
            raise Exception("poly_data input is not a dictionary")
        
        #now on to extracting the data of interest
        print('Assigning element attribute IDs...')
        for i, key in enumerate(poly_data):
            poly_x=poly_data[key][0]#polygon x coordinates
            poly_y=poly_data[key][1]#polygon y coordinates
            inside = isinpolygon(np.array(elm_xyz[0]),
                                 np.array(elm_xyz[2]),
                                 (poly_x,poly_y))
            material_no[inside]=i+1
                            
        self.zone = material_no
        return material_no
    
    def assign_zone_3D(self,volume_data):
        """Assign material/region assocations with certain elements in the mesh 
        say if you have an area you'd like to forward model. 
        ***3D ONLY***
        
        Parameters
        -----------
        volume_data : dict
            Each key contains columns of polygon data for each volume in 
            the form (polyx, polyy, polyz), the polygon data should be the 
            face coordinates which bound the volume.
                        
        Returns
        -----------
        material_no : numpy.array
            Element associations starting at 1. So 1 for the first region 
            defined in the region_data variable, 2 for the second region 
            defined and so on. If the element can't be assigned to a region
            then it'll be left at 0. 
        """
        no_elms=self.num_elms#number of elements 
        elm_xy=self.elm_centre#centriods of mesh elements 
        material_no=np.zeros(no_elms,dtype=int)#attribute number
        if not isinstance(volume_data,dict):
            raise Exception("poly_data input is not a dictionary")  
        
        for i, key in enumerate(volume_data):
            faces_list = volume_data[key]
            
            inside = isinvolume(elm_xy[0],
                                elm_xy[1],
                                elm_xy[2],
                                faces_list)
            material_no[inside]=i+1
                            
        self.zone = material_no
        return material_no        
        
    def assign_zone_attribute(self,material_no,attr_list,new_key):
        """
        Asssigns values to the mesh which depend on region / material only. E.G 
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
        """
        Applies a function to a mesh by zone number and mesh parameter.
        
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
        
    def reciprocal(self,attr="Resistivity",new_key="Conductivity"):
        """
        Compute reciprocal for a given attribute (ie 1 over that number)
        Example is Resistivity to conductivty conversion. 
        
        Parameters
        ----------
        attr: string
            mesh attribute to compute reciprocal of
        new_key: string
            name of new attribute, 
        """
        self.attr_cache[new_key] = 1/np.array(self.attr_cache[attr])
        self.no_attributes += 1
        
    def computeElmDepth(self,datum_x,datum_y,datum_z, method='bilinear'):
        """ Compute the depth of elements given a datum (or surface).
        
        Parameters
        ----------
        datum_x: array like
            X coordinates of datum 
        datum_y: array like
            Y coordinates of datum, if using 2D mesh then set to 'None' 
        datum_z: array like
            Elevation of datum
        method: str, optional
            Method of interpolation used to compute cell depths for a 3D mesh.
            Ignored for 2D meshes. 
            - 'bilinear' - (default) binlinear interpolation
            - 'idw' - inverse distance wieghting
            - 'nearest' - nearest neighbour interpolation
        """
        #formalities and error checking
        if datum_y is None: # set up y column if not in use
            datum_y = [0]*len(datum_x)
        if len(datum_x) != len(datum_y) and len(datum_x) != len(datum_z):
            raise ValueError("Mis match in array dimensions for datum x y z coordinates.")
        datum_x = np.array(datum_x)
        datum_y = np.array(datum_y)
        datum_z = np.array(datum_z)
        if self.ndims == 2: # use 2D interpolation
            elm_x = np.array(self.elm_centre[0])
            elm_z = np.array(self.elm_centre[2])
            min_idx = np.argmin(datum_x)
            max_idx = np.argmax(datum_x)
            Z = np.interp(elm_x,datum_x,datum_z,left=datum_y[min_idx],right=datum_y[max_idx])
            depth = Z - elm_z
            self.attr_cache['depths'] = depth
            self.no_attributes += 1
            return depth
        if self.ndims == 3: # use 3D interpolation
            elm_x = np.array(self.elm_centre[0])
            elm_y = np.array(self.elm_centre[1])
            elm_z = np.array(self.elm_centre[2])  
            #use interpolation to work out depth to datum 
            if method == 'bilinear':
                Z = interp.bilinear(elm_x, elm_y, datum_x, datum_y, datum_z)
            elif method == 'idw':
                Z = interp.idw(elm_x, elm_y, datum_x, datum_y, datum_z)
            elif method == 'nearest':
                Z = interp.nearest(elm_x, elm_y, datum_x, datum_y, datum_z)
            else:
                avail_methods = ['bilinear','idw','nearest']
                raise NameError("Unknown interpolation method, available methods are %s"%str(avail_methods))
            depth = Z - elm_z
            self.attr_cache['depths'] = depth # add cell depths to attribute cache
            self.no_attributes += 1
            return depth
            
    def move_elec_nodes(self, new_x, new_y, new_z, debug=True):
        """
        Move the electrodes to different nodes which are close to the given coordinates. 
        This is useful for timelapse surveys where the electrodes move through time, 
        ideally this is implimented on a mesh which is refined near the surface. If 
        no nodes are assigned to mesh, a mesh.e_nodes variable is created.  
        
        Parameters
        ------------
        new_x : array like
            new electrode x coordinates 
        new_y : array like
            new electrode y coordinates, if 2D mesh let new_y=None and the array will automatically be
            assigned an array of zeros. 
        new_z : array-like
            new electrode z coordinates 
        debug : bool, optional
            Controls if any changes to electrode nodes will be output to console. 
        Returns
        ------------
        node_in_mesh : np array
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
        
    def write_dat(self,file_path='mesh.dat', param=None, zone=None):
        """
        Write a mesh.dat kind of file for mesh input for R2. R2 takes a mesh
        input file for triangle meshes, so this function is only relevant for
        triangle meshes.
        
        Parameters
        ------------
        file_path : string, optional
            Path to the file. By default 'mesh.dat' is saved in the working directory. 
        zone : array like, optional
            An array of integers which are assocaited with regions/materials in the mesh. 
            Useful in the case of an inversion which has a boundary constraint. 
            You can use assign_zone to give a zone/material number to the mesh, and pass that 
            as the zone argument. 
        param : array-like, optional
            Array of parameter number. Set a parameter number to zero fixed its
            conductivity to the starting conductivity.
        
        Notes
        ----------
        mesh.dat like file written to file path. 
        ***IMPORTANT***
        R2/FORTRAN indexing starts at one, in python indexing natively starts at 0
        so when writing mesh.dat we need to check that the node indexes match up correctly.
        """
        if not isinstance(file_path,str):
            raise TypeError("expected string argument for file_path")
        ### write data to mesh.dat kind of file ###
        #open mesh.dat for input      
        fid=open(file_path, 'w')
        
        #write to mesh.dat total num of elements and nodes
        if self.ndims==3:
            fid.write('%i %i 1 0 4\n'%(self.num_elms,self.num_nodes))
        else:
            fid.write('%i %i\n'%(self.num_elms,self.num_nodes))

        #compute zones if present 
        if zone  is None:
            zone = np.ones(self.num_elms, dtype=int) # default zone = 1 
        else:
            if len(zone) != self.num_elms:
                raise IndexError("the number of zone parameters does not match the number of elements")
            elif min(zone) == 0:
                zone = np.array(zone,dtype=int)+1 # as fortran indexing starts at 1, not 0 we must add one to the array if min ==0 
        
        if param  is None:
            param = 1 + np.arange(self.num_elms) # default one parameter per element
        else:
            if len(param) != self.num_elms:
                raise IndexError("the number of parameters does not match the number of elements")
        
        #write out elements         
        no_verts = self.Type2VertsNo()
        for i in range(self.num_elms):
            elm_no=i+1
            fid.write("%i "%elm_no)
            [fid.write("%i "%(self.con_matrix[k][i]+1)) for k in range(no_verts)]
            fid.write("%i %i\n"%(param[i],zone[i]))
    
        #now add nodes
        x_coord = self.node_x
        y_coord = self.node_y
        z_coord = self.node_z
        if self.ndims==3:
            for i in range(self.num_nodes):
                ni_no=i+1
                fid.write("%i %6.3f %6.3f %6.3f\n"%#node number, x coordinate, y coordinate, z coordinate
                          (ni_no,
                           x_coord[i],
                           y_coord[i],
                           z_coord[i]))
            fid.write('1')
        else:
            for i in range(self.num_nodes):
                ni_no=i+1
                fid.write("%i %6.3f %6.3f\n"%#node number, x coordinate, y coordinate
                          (ni_no,
                           x_coord[i],
                           z_coord[i]))

        fid.close()#close the file 
        print('written mesh.dat file to \n%s'%file_path)

    def write_vtk(self,file_path="mesh.vtk", title=None):
        """
        Writes a vtk file for the mesh object, everything in the attr_cache
        will be written to file as attributes. We suggest using Paraview 
        to display the mesh outside of pyR2. It's fast and open source :). 
        
        Parameters
        ------------
        file_path : string, optional
            Maps where python will write the file, if left as `default` then mesh.vtk
            will be written the current working directory. 
        title : string, optional
            Header string written at the top of the vtk file .
        """
        #open file and write header information    
        fh = open(file_path,'w')
        fh.write("# vtk DataFile Version 3.0\n")
        if title == None:
            try:
                title = self.mesh_title
            except AttributeError:
                title = "output from R2 gui meshTools module"
        fh.write(title+"\n")
        fh.write("ASCII\nDATASET UNSTRUCTURED_GRID\n")
        #define node coordinates
        fh.write("POINTS %i double\n"%self.num_nodes)
        for i in range(self.num_nodes):
            if self.ndims == 2:
                fh.write("%8.6f\t%8.6f\t%8.6f\n"%(self.node_x[i],self.node_z[i],self.node_y[i])) # there is no Z in 2D
            elif self.ndims == 3:
                fh.write("%8.6f\t%8.6f\t%8.6f\n"%(self.node_x[i],self.node_y[i],self.node_z[i]))
        #define the connection matrix    
        no_verts = self.Type2VertsNo()
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
        for i,key in enumerate(self.attr_cache):
            fh.write("SCALARS %s double 1\n"%key)
            fh.write("LOOKUP_TABLE default\n")
            [fh.write("%8.6f "%self.attr_cache[key][j]) for j in range(self.num_elms)]
            fh.write("\n")
        
        #finish writing
        fh.write("POINT_DATA %i"%self.num_nodes)        
        fh.close()
    

    def write_attr(self,attr_key,file_name='_res.dat',file_path='default'):
        """ 
        Writes a attribute to a _res.dat type file. file_name entered
        seperately because it will be needed for the R2 config file.
        The reason for this function is so you can write a forward model 
        parameter input file. 
        
        Parameters
        ----------
        attr_key : string
            Key identifying the attr to be written in the mesh object attr_cache.
        file_name : string, optional
            Name of the _res.dat type file.
        file_path : string, optional
            Directory to which the file will be saved in, if left as none then the
            file will be written in the current working directory.
        """
        #formality checks 
        if len(file_name)>15:
            raise NameError("File name for _res.dat type file cannot be longer than 15 characters")
            
        if isinstance(file_name,str)==False or isinstance(file_path,str) == False:
            raise NameError("file_name and file_path arguments must be strings")
        
        if file_path == 'default':#no directory given then ignore file path input
            file_path = file_name
        else:#reassign file_path to full path including the name
            file_path = os.path.join(file_path,file_name)
        
        #the format of the _res.dat file is such that
        #| x coordinate | y coordinate | value | log(value) | 
        fh = open(file_path,'w')#open file handle 
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
        path = drive_letter+'\Program Files'
        contents = os.listdir(path)
        found = False
        for i,pname in enumerate(contents):
            if pname.find("ParaView") != -1:
                para_dir = os.path.join(path,pname)
                found = True
                break
    
        if not found:#try looking in x86 porgram files instead
            path = drive_letter+'\Program Files (x86)'
            contents = os.listdir(path)
            for i,pname in enumerate(contents):
                if pname.find("ParaView") != -1:
                    para_dir = os.path.join(path,pname)
                    found = True
                    break
        
        if not found:
            return False, 'n/a' 
        else:
            return True, os.path.join(para_dir,'bin\paraview.exe')
        #the string output can be run in the console if it is enclosed in speech
        #marks , ie <"C/program files/ParaView5.X/bin/paraview.exe">
        
    def paraview(self,fname='TRIP4Dmesh.vtk',loc=None):
        """
        Show mesh in paraview 
        
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
        op_sys = platform.system()#find kernal type
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
                os.popen(cmd_line+' '+fname)
            except PermissionError:
                print("Windows has blocked launching paraview, program try running as admin")      
        else:
            Popen(['paraview', fname])
            
    def quadMeshNp(self, topo=None):
        """Convert mesh nodes into x column indexes in the case of quad meshes. 
        Does not currently support changes in electrode elevation! 
        
        Returns
        ----------
        colx: list
            X column indexes for quad mesh 
        """
        if int(self.cell_type[0])==8 or int(self.cell_type[0])==9:#elements are quads
            pass
        else:
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
        
#%% triangle centriod 
def tri_cent(p,q,r):
    """
    Compute the centre coordinates for a 2d triangle given the x,y coordinates 
    of the vertices.
            
    Parameters
    ----------
    p : tuple,list,np array
        Coordinates of triangle vertices in the form (x,y).
    q : tuple,list,np array
        Coordinates of triangle vertices in the form (x,y).
    r : tuple,list,np array
        Coordinates of triangle vertices in the form (x,y).
            
    Returns
    ----------
    coordinates : tuple
        In the format (x,y).    
    """
    Xm=(p[0]+q[0])/2
    Ym=(p[1]+q[1])/2
    k=2/3
    Xc=r[0]+(k*(Xm-r[0]))
    Yc=r[1]+(k*(Ym-r[1]))
    return(Xc,Yc)
    
#%% import a vtk file 
def vtk_import(file_path='mesh.vtk',parameter_title='default'):
    """
    Imports a mesh file into the python workspace, can have triangular, quad or tetraheral shaped elements.
            
    Parameters
    ----------
    file_path : string, optional
        File path to mesh file. Note that a error will occur if the file format is not as expected.
    parameter_title : string, optional
        Name of the parameter table in the vtk file, if left as default the first look up table found will be returned 
        also note that all parameters will be imported. Just the title highlights which one the mesh object will use as 
        default cell attribute. 
            
    Returns
    -------
    mesh : class 
        a <pyR2> mesh class 
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
        no_nodes=int(node_info[1])
    except IndexError:#if we get this then there is a white space between the node info and header lines
        node_info=fid.readline().strip().split()#read line 5
        no_nodes=int(node_info[1])      
    if no_nodes == 0: 
        raise ImportError("No nodes in vtk file to import! Aborting... ")
    #now read in node data
    x_coord=[]#make lists for each of the relevant parameters for each node
    y_coord=[]
    z_coord=[]
    node_num=[]
    for i in range(no_nodes):
        try:
            coord_data=fid.readline().strip().split()
            x_coord.append(float(coord_data[0]))
            y_coord.append(float(coord_data[1]))
            z_coord.append(float(coord_data[2]))
        except:# ValueError:
            coord_data=fid.readline()
            x_coord.append(float(coord_data[0:12])) # retrive fixed width columns if cannot parse as split strings
            y_coord.append(float(coord_data[12:24]))
            z_coord.append(float(coord_data[24:36]))
        node_num.append(i)
    
    #now read in element data
    #print("importing mesh element info...")
    elm_info=fid.readline().strip().split()#read line with cell data
    try:
        no_elms=int(elm_info[1])
    except IndexError: # quick bug fix
        elm_info=fid.readline().strip().split()#read line with cell data
        no_elms=int(elm_info[1])
    
    if no_elms ==0: 
        raise ImportError("No elements in vtk file to import!")
    
    no_pts=[]#assign lists to nodes 
    node1=[]
    node2=[]
    node3=[]
    node4=[]
    node5=[]
    node6=[]
    node7=[]
    node8=[]
    #node9=[]
    elm_num=[]
    centriod_x=[]#list will contain the centre points of elements 
    centriod_y=[]
    centriod_z=[]
    areas=[]#areas of cells (might be useful in the future)
    ignored_cells=0
    #import element data ... expects triangles or quads 
    for i in range(no_elms):
        elm_data=fid.readline().strip().split()
        if int(elm_data[0])==3:
            if i==0:
                #print("triangular elements detected")
                vert_no=3
            no_pts.append(int(elm_data[0]))
            #nodes
            node1.append(int(elm_data[1]))
            node2.append(int(elm_data[2]))
            node3.append(int(elm_data[3]))
            elm_num.append(i+1)
            #find the centriod of the element for triangles
            n1=(x_coord[int(elm_data[1])],y_coord[int(elm_data[1])])#in vtk files the 1st element id is 0 
            n2=(x_coord[int(elm_data[2])],y_coord[int(elm_data[2])])
            n3=(x_coord[int(elm_data[3])],y_coord[int(elm_data[3])])
            xy_tuple=tri_cent(n1,n2,n3)#actual calculation
            centriod_x.append(xy_tuple[0])
            centriod_y.append(xy_tuple[1])
            #find area of element (for a triangle this is 0.5*base*height)
            base=(((n1[0]-n2[0])**2) + ((n1[1]-n2[1])**2))**0.5
            mid_pt=((n1[0]+n2[0])/2,(n1[1]+n2[1])/2)
            height=(((mid_pt[0]-n3[0])**2) + ((mid_pt[1]-n3[1])**2))**0.5
            areas.append(0.5*base*height)
        elif int(elm_data[0])==4:
            if i==0:
                vert_no=4
            no_pts.append(int(elm_data[0]))
            #nodes
            node1.append(int(elm_data[1]))
            node2.append(int(elm_data[2]))
            node3.append(int(elm_data[3]))
            node4.append(int(elm_data[4]))
            elm_num.append(i+1)
            #assuming element centres are the average of the x - y coordinates for the quad
            n1=(x_coord[int(elm_data[1])],y_coord[int(elm_data[1])],z_coord[int(elm_data[1])])#in vtk files the 1st element id is 0 
            n2=(x_coord[int(elm_data[2])],y_coord[int(elm_data[2])],z_coord[int(elm_data[2])])
            n3=(x_coord[int(elm_data[3])],y_coord[int(elm_data[3])],z_coord[int(elm_data[3])])
            n4=(x_coord[int(elm_data[4])],y_coord[int(elm_data[4])],z_coord[int(elm_data[4])])
            centriod_x.append(np.mean((n1[0],n2[0],n3[0],n4[0])))
            centriod_y.append(np.mean((n1[1],n2[1],n3[1],n4[1])))
            centriod_z.append(np.mean((n1[2],n2[2],n3[2],n4[2])))
            #finding element areas, base times height.  
            elm_len=abs(n2[0]-n1[0])#element length
            elm_hgt=abs(n2[1]-n3[1])#element hieght
            areas.append(elm_len*elm_hgt)
        elif int(elm_data[0])==8: # this following code is getting silly in how long it is. Need to work on a more efficent way
            if i==0:
                vert_no=8
            no_pts.append(int(elm_data[0]))
            #nodes
            node1.append(int(elm_data[1]))
            node2.append(int(elm_data[2]))
            node3.append(int(elm_data[3]))
            node4.append(int(elm_data[4]))
            node5.append(int(elm_data[5]))
            node6.append(int(elm_data[6]))
            node7.append(int(elm_data[7]))
            node8.append(int(elm_data[8]))
            #assuming element centres are the average of the x - y coordinates for the quad
            n1=(x_coord[int(elm_data[1])],y_coord[int(elm_data[1])],z_coord[int(elm_data[1])])#in vtk files the 1st element id is 0 
            n2=(x_coord[int(elm_data[2])],y_coord[int(elm_data[2])],z_coord[int(elm_data[2])])
            n3=(x_coord[int(elm_data[3])],y_coord[int(elm_data[3])],z_coord[int(elm_data[3])])
            n4=(x_coord[int(elm_data[4])],y_coord[int(elm_data[4])],z_coord[int(elm_data[4])])
            n5=(x_coord[int(elm_data[5])],y_coord[int(elm_data[5])],z_coord[int(elm_data[5])]) 
            n6=(x_coord[int(elm_data[6])],y_coord[int(elm_data[6])],z_coord[int(elm_data[6])])
            n7=(x_coord[int(elm_data[7])],y_coord[int(elm_data[7])],z_coord[int(elm_data[7])])
            n8=(x_coord[int(elm_data[8])],y_coord[int(elm_data[8])],z_coord[int(elm_data[8])])
            centriod_x.append(np.mean((n1[0],n2[0],n3[0],n4[0],n5[0],n6[0],n7[0],n8[0])))
            centriod_y.append(np.mean((n1[1],n2[1],n3[1],n4[1],n5[1],n6[1],n7[1],n8[1])))
            centriod_z.append(np.mean((n1[2],n2[2],n3[2],n4[2],n5[2],n6[2],n7[2],n8[2])))
            #estimate element VOLUMES, base area times height.  
            elm_len=abs(n2[0]-n1[0])#element length
            elm_width = abs(n1[1]-n3[1])
            elm_thick=abs(n5[2]-n1[2])#element hieght
            areas.append(elm_len*elm_width*elm_thick)
            
        else: 
            warnings.warn("WARNING: unkown cell type encountered!")
            ignored_cells+=1
    #compile some information   
    
    if len(centriod_z)==0:#check if mesh is 2D 
        centriod_z==[0]*len(centriod_x)
    if sum(z_coord)==0:#then mesh is 2D and node y and node z, centriod y and centriod z columns need swapping so they work in mesh tools 
        temp_y = y_coord
        y_coord = z_coord
        z_coord = temp_y
        temp_y = centriod_y
        centriod_y = centriod_z
        centriod_z = temp_y
    
    centriod=(centriod_x,centriod_y,centriod_z)#centres of each element in form (x...,y...,z...)
    if vert_no==3:
        node_maps=(node1,node2,node3)
    elif vert_no==4:
        node_maps=(node1,node2,node3,node4)  
    elif vert_no==8:
        node_maps=(node1,node2,node3,node4,node5,node6,node7,node8)
        
    if ignored_cells>0:
        print("%i cells ignored in the vtk file"%ignored_cells)
    
    cell_attr_dump=fid.readlines()#reads the last portion of the file
    #finished reading the file
    
    #find cell types
    for i,line_info in enumerate(cell_attr_dump):
        if line_info.find("CELL_TYPES") == 0:
            cell_type = [int(k) for k in cell_attr_dump[i+1].strip().split()]
            break
    
    fid.close()
    #print("reading cell attributes...")
    # read through cell attributes to find the relevant parameter table?
    
    #find scalar values in the vtk file
    num_attr = 0
    attr_dict = {}
    #found = False # boolian if we have found the parameter of interest
    for i,line_info in enumerate(cell_attr_dump):
        if line_info.find("SCALARS") == 0:
            attr_title = line_info.split()[1]
            #check look up table
            if cell_attr_dump[i+1].split()[1] != "default":
                warnings.warn("unrecognised lookup table type")
            values=[float(k) for k in cell_attr_dump[i+2].split()]
            attr_dict[attr_title] = values
            if num_attr ==0:# primary attribute defaults to the first attribute found
                parameter_title = attr_title
                values_oi = values
            if attr_title == parameter_title:#then its the parameter of interest that the user was trying extract
                #found = True
                values_oi = values        
            num_attr += 1
    
    #put in fail safe if no attributes are found        
    if num_attr == 0:
        print("no cell attributes found in vtk file")
        attr_dict = {"no attributes":[float("nan")]*no_elms}
        values_oi= [0]*no_elms
        parameter_title = "n/a"
    #print("finished importing mesh.\n")
    #information in a dictionary, this is easier to debug than an object in spyder: 
    mesh_dict = {'num_nodes':no_nodes,#number of nodes
            'num_elms':no_elms,#number of elements 
            'node_x':x_coord,#x coordinates of nodes 
            'node_y':y_coord,#y coordinates of nodes
            'node_z':z_coord,#z coordinates of nodes 
            'node_id':node_num,#node id number 
            'elm_id':elm_num,#element id number 
            'num_elm_nodes':no_pts,#number of points which make an element
            'node_data':node_maps,#nodes of element vertices
            'elm_centre':centriod,#centre of elements (x,y)
            'elm_area':areas,#area of each element (or volume)
            'cell_type':cell_type,
            'parameters':values_oi,#the values of the attributes given to each cell 
            'parameter_title':parameter_title,
            'cell_attributes':attr_dict,
            'dict_type':'mesh_info',
            'original_file_path':file_path} 
    mesh = Mesh.mesh_dict2class(mesh_dict)#convert to mesh object
#    print(mesh.attr_cache.keys())
    try:
        if mesh.ndims==2:
            mesh.add_sensitivity(mesh.attr_cache['Sensitivity(log10)'])
        else:
            mesh.add_sensitivity(mesh.attr_cache['Sensitivity_map(log10)'])
    except:
        print('no sensitivity')

    mesh.mesh_title = title
    return mesh

#%% import mesh from native .dat format
def dat_import(file_path='mesh.dat'):
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
    node_map = np.array([[0]*numel]*npere,dtype=int)
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
    #iniate mesh class 
    #compute each cell area and centriod 
    areas = [0]*numel
    centriod_x = [0]*numel
    centriod_y = [0]*numel
    centriod_z = [0]*numel
    for i in range(numel):
        if npere==3:
            n1=(node_x[node_map[0][i]],node_z[node_map[0][i]])#in vtk files the 1st element id is 0 
            n2=(node_x[node_map[1][i]],node_z[node_map[1][i]])
            n3=(node_x[node_map[2][i]],node_z[node_map[2][i]])
            xy_tuple=tri_cent(n1,n2,n3)#actual calculation
            centriod_x[i] = xy_tuple[0]
            centriod_z[i] = xy_tuple[1]
            #find area of element (for a triangle this is 0.5*base*height)
            base=(((n1[0]-n2[0])**2) + ((n1[1]-n2[1])**2))**0.5
            mid_pt=((n1[0]+n2[0])/2,(n1[1]+n2[1])/2)
            height=(((mid_pt[0]-n3[0])**2) + ((mid_pt[1]-n3[1])**2))**0.5
            areas[i] = 0.5*base*height
        else:
            x_vec = [node_x[node_map[j][i]] for j in range(npere)]
            y_vec = [node_y[node_map[j][i]] for j in range(npere)]
            z_vec = [node_z[node_map[j][i]] for j in range(npere)]
            centriod_x[i] = sum(x_vec)/npere
            centriod_y[i] = sum(y_vec)/npere
            centriod_z[i] = sum(z_vec)/npere
            #dont compute area as it is not needed 
    
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
        
    mesh = Mesh(num_nodes = numnp,#number of nodes
                 num_elms = numel,#number of elements 
                 node_x = node_x,#x coordinates of nodes 
                 node_y = node_y,#y coordinates of nodes
                 node_z = node_z,#z coordinates of nodes 
                 node_id= node_id,#node id number 
                 elm_id=elm_no,#element id number 
                 node_data=node_map,#nodes of element vertices
                 elm_centre= (centriod_x,centriod_y,centriod_z),#centre of elements (x,y)
                 elm_area = areas,#area of each element
                 cell_type = [cell_type],#according to vtk format
                 cell_attributes = zone,#the values of the attributes given to each cell, we dont have any yet 
                 atribute_title='zone')#what is the attribute? we may use conductivity instead of resistivity for example
    
    mesh.add_attribute(zone,'zone')
    
    return mesh 
        
           
#%% Read in E4D / tetgen mesh
def tetgen_import(file_path):
    """Import Tetgen mesh into <pyR2>. This isa little different from other 
    imports as the mesh is described by several files. From pyR2's perspective
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
    
    #next we need the connection matrix to describe each of the tetrahedral e
    #elements
    file_path2 = file_path.replace('.node','.ele')
    fh = open(file_path2,'r')# read in element file  
    header = fh.readline() # read in header line 
    numel = int(header.split()[0]) # number of elements 
    npere = 4 # number of nodes per element, in this case E4D uses tetrahedral 
    #meshes so its always going to be 4. 
    
    
    node_map = np.array([[0]*numel]*npere,dtype=int) # connection matrix mapping elements onto nodes 
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
    
    #calculate element centres       
    areas = [0]*numel
    centriod_x = [0]*numel
    centriod_y = [0]*numel
    centriod_z = [0]*numel    
    for i in range(numel):
        x_vec = [node_x[node_map[j][i]] for j in range(npere)]
        y_vec = [node_y[node_map[j][i]] for j in range(npere)]
        z_vec = [node_z[node_map[j][i]] for j in range(npere)]
        centriod_x[i] = sum(x_vec)/npere
        centriod_y[i] = sum(y_vec)/npere
        centriod_z[i] = sum(z_vec)/npere
            
    
    #create mesh instance 
    mesh = Mesh(num_nodes = numnp,#number of nodes
             num_elms = numel,#number of elements 
             node_x = node_x,#x coordinates of nodes 
             node_y = node_y,#y coordinates of nodes
             node_z = node_z,#z coordinates of nodes 
             node_id= node_id,#node id number 
             elm_id=elm_no,#element id number 
             node_data=node_map,#nodes of element vertices
             elm_centre= (centriod_x,centriod_y,centriod_z),#centre of elements (x,y)
             elm_area = areas,#area of each element
             cell_type = [10],#according to vtk format
             cell_attributes = zone,#the values of the attributes given to each cell, we dont have any yet 
             atribute_title='zone')#what is the attribute? 
    
    mesh.add_attribute(zone,'zone')
    
    return mesh
        
#%% build a quad mesh        
def quad_mesh(elec_x, elec_z, elec_type = None, elemx=4, xgf=1.5, yf=1.1, ygf=1.25, doi=-1, pad=2, 
              surface_x=None,surface_z=None):
    """Creates a quaderlateral mesh given the electrode x and y positions. Function
    relies heavily on the numpy package.
            
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
    yf : float, optional
         Y factor multiplier in the fine zone.
    ygf : float, optional
         Y factor multiplier in the coarse zone.
    doi : float (m), optional 
         Depth of investigation (if left as -1 = half survey width).
    pad : int, optional
         X padding outside the fine area (tipicaly twice the number of elements between electrodes).
    surface_x: array like, optional
        Default is None. x coordinates of extra surface topography points for the generation of topography in the quad mesh
    surface_y: array like, optional
        Default is None. y coordinates of extra surface topography points for the generation of topography in the quad mesh. Note
        an error will be returned if len(surface_x) != len(surface_y)
        
            
    Returns
    -------
    Mesh : class
        Mesh object 
    meshx : numpy.array
        Mesh x locations for R2in file.
    meshy : numpy.array
        Mesh y locations for R2in file (ie node depths).
    topo : numpy.array
        Topography for R2in file.
    elec_node : numpy.array
        x columns where the electrodes are. 
    """
    #formalities, error check
    if elemx < 4:
        print('elemx too small, set up to 4 at least')
        elemx = 4
        
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
            Ex=np.array(elec_x)[surface_idx]
            Ey=np.array(elec_z)[surface_idx]
        elif len(surface_idx)== 0 and len(surface_x)>0:
            #case where you have surface topography but no surface electrodes 
            Ex=np.array(surface_x)
            Ey=np.array(surface_z)
            elec=np.c_[Ex,Ey]
        elif len(surface_idx)== 0:
            #fail safe if no surface electrodes are present to generate surface topography 
            Ex=np.array([elec_x[np.argmin(elec_x)],elec_x[np.argmax(elec_x)]])
            Ey=np.array([elec_z[np.argmax(elec_z)],elec_z[np.argmax(elec_z)]])
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
            
    # create meshy
    if doi == -1:
        if bh_flag:
            doi = abs(min(elec_z))
        else:
            doi = np.abs(elec[0,0]-elec[-1,0])/2
        
#    dyy = espacing/(elemx*4)
    meshy = [0]
    dyy = 0.05
    for i in range(100):
        meshy.append(meshy[-1]+dyy*yf)
        dyy = dyy*yf
        if meshy[-1] > doi:
            break
    elemy = len(meshy)
    elemy2 = int(elemy/2)
    yy = np.ones(elemy2)*meshy[-1]
    for i in range(1, elemy2):
        yy[i] = yy[i-1]+dyy*ygf
        dyy = dyy*ygf
    meshy = np.r_[meshy, yy[1:]]
    
    #insert borehole electrodes? if we have boreholes / buried electrodes 
    if bh_flag:
        meshx = np.unique(np.append(meshx,bh[:,0]))
        
    # create topo
    if bh_flag: # only use surface electrodes to make the topography if buried electrodes present
        X = np.append(Ex,surface_x) 
        Y = np.append(Ey,surface_z)
        idx = np.argsort(X)
        topo = np.interp(meshx, X[idx], Y[idx])
    else: # all electrodes are assumed to be on the surface 
        X = np.append(elec[:,0],surface_x)
        Y = np.append(elec[:,1],surface_z)
        idx = np.argsort(X)
        topo = np.interp(meshx, X[idx], Y[idx])
    
    if bh_flag:
        #insert y values of boreholes, normalised to topography
        norm_bhy = np.interp(bh[:,0], elec[:,0], elec[:,1]) - bh[:,1]
        meshy = np.unique(np.append(meshy,norm_bhy))
    
    ###
    #find the columns relating to the electrode nodes? 
    temp_x = meshx.tolist()
    temp_y = meshy.tolist()
    elec_node_x=[temp_x.index(elec_x[i])+1 for i in range(len(elec_x))]#add 1 becuase of indexing in R2. 
    
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
            raise Exception("There was a problem indexing the meshy values for electrode number %i"%i)
                
    elec_node = [elec_node_x,elec_node_y.T.tolist()[0]]
    
    #print some warnings for debugging 
    if len(topo)!=len(meshx):
        warnings.warn("Topography vector and x coordinate arrays not the same length! ")
    elif len(elec_node_x)!=len(elec_x):
        warnings.warn("Electrode node vector and number of electrodes mismatch! ")
     
    # what is the number of regions? (elements)
    no_elms=(len(meshx)-1)*(len(meshy)-1)
    no_nodes=len(meshx)*len(meshy)
    
    # compute node mappins (connection matrix)
    y_dim=len(meshy)
    fnl_node=no_nodes-1
    
    node_mappins=(np.arange(0,fnl_node-y_dim),
                  np.arange(y_dim,fnl_node),
                  np.arange(y_dim+1,fnl_node+1),
                  np.arange(1,fnl_node-y_dim+1))
    
    del_idx=np.arange(y_dim-1,len(node_mappins[0]),y_dim)#the above has too many indexes at the changeover of columns so some need deleting
    
    node_mappins = [list(np.delete(node_mappins[i],del_idx)) for i in range(4)]#delete excess node placements
    #compute node x and y  (and z)
    node_x,node_z=np.meshgrid(meshx,meshy)
    #account for topography in the y direction 
    node_z = [topo-node_z[i,:] for i in range(y_dim)]#list comprehension to add topography to the mesh, (could use numpy here??)
    node_z=np.array(node_z).flatten(order='F')
    node_x=node_x.flatten(order='F')
    node_y=np.array([0]*len(node_x))
    
    #compute element centres and areas
    centriod_x=[]
    centriod_z=[]
    areas=[]
    for i in range(no_elms):
        #assuming element centres are the average of the x - y coordinates for the quad
        n1=(node_x[int(node_mappins[0][i])],node_z[int(node_mappins[0][i])])#in vtk files the 1st element id is 0 
        n2=(node_x[int(node_mappins[1][i])],node_z[int(node_mappins[1][i])])
        n3=(node_x[int(node_mappins[2][i])],node_z[int(node_mappins[2][i])])
        n4=(node_x[int(node_mappins[3][i])],node_z[int(node_mappins[3][i])])
        centriod_x.append(np.mean((n1[0],n2[0],n3[0],n4[0])))
        centriod_z.append(np.mean((n1[1],n2[1],n3[1],n4[1])))
        #finding element areas, base times height.  
        elm_len=abs(n2[0]-n1[0])#element length
        elm_hgt=abs(n2[1]-n3[1])#element hieght
        areas.append(elm_len*elm_hgt)
    centriod_y = [0]*len(centriod_x)
    
    #make mesh class    
    mesh = Mesh(no_nodes,
                    no_elms,
                    node_x,
                    node_y,
                    node_z,
                    list(np.arange(0,no_nodes)),
                    list(np.arange(0,no_elms)),
                    node_mappins,
                    (centriod_x,centriod_y,centriod_z),
                    areas,
                    [9],
                    [0]*no_elms,
                    'no attribute')
    
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

    return mesh,meshx,meshy,topo,elec_node


#%% build a triangle mesh - using the gmsh wrapper
def tri_mesh(elec_x, elec_z, elec_type=None, geom_input=None,keep_files=True, 
             show_output=True, path='exe', dump=print, whole_space=False, **kwargs):
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
        <pyR2> mesh class
        
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
    
    if not os.path.isfile(os.path.join(ewd,'gmsh.exe')):
        raise Exception("No gmsh.exe exists in the exe directory!")
    
    #make .geo file
    file_name="mesh"
    if not whole_space:#by default create survey with topography 
        node_pos = gw.genGeoFile([elec_x,elec_z], elec_type, geom_input,
                             file_path=file_name,**kwargs)
    elif whole_space:
        print("Whole space problem")
        node_pos = gw.gen_2d_whole_space([elec_x,elec_z], geom_input = geom_input, 
                                         file_path=file_name)    
    
    # handling gmsh
    if platform.system() == "Windows":#command line input will vary slighty by system 
        cmd_line = ewd+'\gmsh.exe '+file_name+'.geo -2'
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
        cmd_line = ['wine',ewd+'/gmsh.exe', file_name+'.geo', '-2']
        
    if show_output: 
        p = Popen(cmd_line, stdout=PIPE, shell=False)#run gmsh with ouput displayed in console
        while p.poll() is None:
            line = p.stdout.readline().rstrip()
            if line.decode('utf-8') != '':
                dump(line.decode('utf-8'))
    else:
        call(cmd_line)#run gmsh 
        
    #convert into mesh.dat 
    mesh_dict = gw.msh_parse(file_path = file_name+'.msh') # read in mesh file
    mesh = Mesh.mesh_dict2class(mesh_dict) # convert output of parser into an object
    #mesh.write_dat(file_path='mesh.dat') # write mesh.dat - disabled as handled higher up in the R2 class 
    
    if keep_files is False: 
        os.remove(file_name+".geo");os.remove(file_name+".msh")

    mesh.add_e_nodes(node_pos-1)#in python indexing starts at 0, in gmsh it starts at 1 
    
    # point at the surface
    xsurf = []
    zsurf = []
    for x, z, t in zip(elec_x, elec_z, elec_type):
        if t == 'electrode': # surface electrode
            xsurf.append(x)
            zsurf.append(z)
    if 'surface' in geom_input.keys():
        xsurf = xsurf + list(geom_input['surface'][0])
        zsurf = zsurf + list(geom_input['surface'][1])
    surfacePoints = np.array([xsurf, zsurf]).T
    isort = np.argsort(xsurf)
    mesh.surface = surfacePoints[isort, :]

    
    return mesh#, mesh_dict['element_ranges']


#%% 3D tetrahedral mesh 
def tetra_mesh(elec_x,elec_y,elec_z=None, elec_type = None, keep_files=True, interp_method = 'bilinear',
               surface_refinement = None, mesh_refinement = None,show_output=True, 
               path='exe', dump=print,whole_space=False, padding=20, search_radius = 10,
               **kwargs):
    """ 
    Generates a tetrahedral mesh for R3t (with topography). returns mesh3d.dat 
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
        type of electrode see Notes in genGeoFile in gmshWrap.py for the format of this list   
    keep_files : boolean, optional
        `True` if the gmsh input and output file is to be stored in the working directory.
    interp_method: string, default ='bilinear' optional
        Interpolation method to translate mesh nodes in the z direction. In other words the method in which topography 
        is appended to the mesh. Here the topography is added to the mesh in post processing. 
        The options are inverse wieghting distance, bilinear interpolation or nearest neighbour. 
        if == 'idw': then provide search_radius.  
    surface_refinement : np.array, optional 
        Numpy array of shape (3,n), should follow the format np.array([x1,x2,x3,...],[y1,y2,y3,...],[z1,z2,z3,...]).
        Allows for extra refinement for the top surface of the mesh. The points are not added to the mesh, but 
        considered in post processing in order to super impose topography on the mesh. 
    mesh_refinement : np.array, optional
        Coming soon ... 
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
    -------
    mesh3d: class
 
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
    avail_methods = ['bilinear','idw','nearest']
    if interp_method not in avail_methods:
        raise NameError("'%s' is an unrecognised interpretation method"%interp_method)
            
    if surface_refinement is not None:
        surf_x = surface_refinement[0,:]
        surf_y = surface_refinement[1,:]
        surf_z = surface_refinement[2,:]
    else:
        surf_x = []
        surf_y = []
        surf_z = []
    
    if elec_type is not None:
        warnings.warn("Borehole electrode meshes still in development!")
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
        for i, key in enumerate(elec_type):
            if key == 'buried':
                bur_elec_x.append(elec_x[i])
                bur_elec_y.append(elec_y[i])
                bur_elec_z.append(elec_z[i])
            if key == 'surface':
                surf_elec_x.append(elec_x[i])
                surf_elec_y.append(elec_y[i])
                surf_elec_z.append(elec_z[i])
        #interpolate in order to normalise buried electrode elevations to 0
        x_interp = np.append(surf_elec_x,surf_x)#parameters to be interpolated with
        y_interp = np.append(surf_elec_y,surf_y)
        z_interp = np.append(surf_elec_z,surf_z)
        
        if interp_method is 'idw': 
            bur_elec_z = np.array(bur_elec_z) - interp.idw(bur_elec_x, bur_elec_y, x_interp, y_interp, z_interp,radius=search_radius)# use inverse distance weighting
        elif interp_method is 'bilinear':
            bur_elec_z = np.array(bur_elec_z) - interp.bilinear(bur_elec_x, bur_elec_y, x_interp, y_interp, z_interp)
        elif interp_method is 'nearest':
            bur_elec_z = np.array(bur_elec_z) - interp.nearest(bur_elec_x, bur_elec_y, x_interp, y_interp, z_interp)
    else:
        surf_elec_x = elec_x 
        surf_elec_y = elec_y 
        if elec_z is None:
            surf_elec_z = np.zeros_like(elec_x)
            elec_z = np.zeros_like(elec_x)
        else:
            surf_elec_z = elec_z 
            elec_z = np.array(elec_z) - np.array(elec_z)#normalise elec_z
         
    #check directories 
    if path == "exe":
        ewd = os.path.join(
                os.path.dirname(os.path.realpath(__file__)),
                path)
        #print(ewd) #ewd - exe working directory 
    else:
        ewd = path # points to the location of the .exe 
        # else its assumed a custom directory has been given to the gmsh.exe 
    
    if not os.path.isfile(os.path.join(ewd,'gmsh.exe')):
        raise Exception("No gmsh.exe exists in the exe directory!")
    
    #make .geo file
    file_name="mesh3d"
    if whole_space:#by default create survey with topography 
        print("Whole space problem")
        raise Exception("Sorry whole space 3D problems are not implimented yet")
        
    else:
        node_pos = gw.box_3d([elec_x,elec_y,elec_z], file_path=file_name, **kwargs)
        
    # handling gmsh
    if platform.system() == "Windows":#command line input will vary slighty by system 
        cmd_line = ewd+'\gmsh.exe '+file_name+'.geo -3'
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
        cmd_line = ['wine',ewd+'/gmsh.exe', file_name+'.geo', '-3']
        
    if show_output: 
        p = Popen(cmd_line, stdout=PIPE, shell=False)#run gmsh with ouput displayed in console
        while p.poll() is None:
            line = p.stdout.readline().rstrip()
            dump(line.decode('utf-8'))
    else:
        call(cmd_line)#run gmsh 
        
    #convert into mesh.dat
    mesh_dict = gw.msh_parse_3d(file_path = file_name+'.msh') # read in 3D mesh file
    mesh = Mesh.mesh_dict2class(mesh_dict) # convert output of parser into an object
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
    if interp_method is 'idw': 
        nodez = interp.idw(node_x, node_y, x_interp, y_interp, z_interp,radius=search_radius)# use inverse distance weighting
    elif interp_method is 'bilinear':# interpolate on a irregular grid, extrapolates the unknown coordinates
        nodez = interp.bilinear(node_x, node_y, x_interp, y_interp, z_interp)
    elif interp_method is 'nearest':
        nodez = interp.nearest(node_x, node_y, x_interp, y_interp, z_interp)
    print('done')
    
    mesh.node_z = np.array(mesh.node_z) + nodez
    node_z = mesh.node_z
    #need to recompute cell centres as well as they will have changed. It's more efficient to recompute them rather than interpolate.
    elm_x = np.array(mesh.elm_centre[0])
    elm_y = np.array(mesh.elm_centre[1])
    elm_z = np.array(mesh.elm_centre[2])
    numel = mesh.num_elms
    npere = 4
    node1 = mesh.con_matrix[0]#indexes for node positions 
    node2 = mesh.con_matrix[1]
    node3 = mesh.con_matrix[2]
    node4 = mesh.con_matrix[3]
    for i in range(numel):
        n1=(node_x[node1[i]],node_y[node1[i]],node_z[node1[i]])#define node coordinates
        n2=(node_x[node2[i]],node_y[node2[i]],node_z[node2[i]])#we have to take 1 off here cos of how python indexes lists and tuples
        n3=(node_x[node3[i]],node_y[node3[i]],node_z[node3[i]])
        n4=(node_x[node4[i]],node_y[node4[i]],node_z[node4[i]])
        elm_x[i] = sum((n1[0],n2[0],n3[0],n4[0]))/npere
        elm_y[i] = sum((n1[1],n2[1],n3[1],n4[1]))/npere
        elm_z[i] = sum((n1[2],n2[2],n3[2],n4[2]))/npere
    
    mesh.elm_centre = (elm_x, elm_y, elm_z)
    #add nodes to mesh
    mesh.add_e_nodes(node_pos-1)#in python indexing starts at 0, in gmsh it starts at 1 
    
    return mesh

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

#%% parser for reading in mesh.dat like file and returning a mesh.     
def dat_import_v1(file_path):
    
    def readMeshDat(fname):
        with open(fname, 'r') as f:
            x = f.readline().split()
        numel = int(x[0])
        elems = np.genfromtxt(fname, skip_header=1, max_rows=numel).astype(int)
        nodes = np.genfromtxt(fname, skip_header=numel+1)
        return elems, nodes
#    
#    def computeCentroid(elems, nodes):
#        elx = nodes[elems,1]
#        ely = nodes[elems,2]
#        cx = np.sum(elx, axis=1)/elems.shape[1]
#        cy = np.sum(ely, axis=1)/elems.shape[1]
#        return np.c_[cx, cy]
    
    def computeCentroid(elems, nodes):
        return np.mean(nodes[elems,:], axis=1)
    
    elems, nodes = readMeshDat(file_path)
    centroids = computeCentroid(elems[:,1:-2]-1, nodes)
    centriod_x, centriod_y = centroids[:,0], centroids[:,1]
    no_nodes = nodes.shape[0]
    no_elms = elems.shape[0]
    x_coord = nodes[:,1]
    z_coord = nodes[:,2]
    node_num = nodes[:,0]
    elm_number = elems[:,0]
    node_data = tuple([list(elems[:,i]-1) for i in range(1, elems.shape[1]-2)])
    if len(node_data) == 3: # triangules
        # using Heron formula
        s = np.sum(elems[:,1:-2], axis=1)/2
        areas = np.sqrt(s*(s-node_data[0])*(s-node_data[1])*(s-node_data[2]))
    elif len(node_data) == 4: # rectuangular
        areas = np.abs(np.array(node_data[0])-np.array(node_data[1]))\
                *np.abs(np.array(node_data[2])-np.array(node_data[3]))
    else: # fill with nan
        areas = np.zeros(elm_number)*np.nan

    mesh = Mesh(num_nodes = no_nodes,#number of nodes
                     num_elms = no_elms,#number of elements 
                     node_x = x_coord,#x coordinates of nodes 
                     node_y = z_coord,#y coordinates of nodes
                     node_z= [0]*no_nodes,#z coordinates of nodes 
                     node_id= node_num,#node id number 
                     elm_id=elm_number,#element id number 
                     node_data=node_data,#nodes of element vertices
                     elm_centre= (centriod_x,centriod_y),#centre of elements (x,y)
                     elm_area = areas,#area of each element
                     cell_type = [5],#according to vtk format
                     cell_attributes = float("nan"),#the values of the attributes given to each cell 
                     atribute_title='none')#what is the attribute? we may use conductivity instead of resistivity for example
    return mesh
#now do mesh.add_e_nodes to add electrode positions to the mesh. 
    
#%% import a custom mesh, you must know the node positions 
def custom_mesh_import(file_path, node_pos=None, flag_3D=False):
    """ 
    Import user defined mesh, currently supports .msh, .vtk and .dat (native to R2/3t)
    format for quad, triangular and tetrahedral meshes. The type of file is guessed from the 
    extension given to the code. 
    
    Parameters
    ---------- 
    file_path: string
        Path to file.
    node_pos: array like, optional
        Array of ints referencing the electrode nodes. If left as none no electrodes 
        will be added to the mesh class. Consider using mesh.move_elec_nodes()
        to add nodes to mesh using their xyz coordinates.
    flag_3D: bool, optional
        Make this true for 3D meshes if importing .msh type. 
        
    Returns
    -------
    mesh: class
        mesh class used in pyR2
        
    """
    if not isinstance(file_path,str):
        raise TypeError("Expected string type argument for 'file_path'")
    
    path,ext = os.path.splitext(file_path)
    if ext == '.vtk':
        mesh = vtk_import(file_path)
    elif ext == '.msh':
        if flag_3D:
            mesh_dict = gw.msh_parse_3d(file_path)
        else:
            mesh_dict = gw.msh_parse(file_path)
        mesh = Mesh.mesh_dict2class(mesh_dict)
    elif ext == '.dat':
        mesh = dat_import(file_path)   
    elif ext == '.node':
        mesh = tetgen_import(file_path)
    else:
        avail_ext = ['.vtk','.msh','.dat','.node / .exe']
        raise ImportError("Unrecognised file extension, available extensions are "+str(avail_ext))
    
    if node_pos is not None:
        mesh.add_e_nodes(np.array(node_pos, dtype=int)) # add electrode nodes to mesh provided by the user
    
    return mesh

#%% ram amount check and is wine installed?. 
#Now for complicated meshes we need alot more RAM. the below function is a os agnostic
#function for returning the amount of total ram. 
# we also need to check wine is installed if running on macOs or linux. 
def systemCheck():
    """
    Performs a simple diagnostic of the system, no input commands needed. System
    info is printed to screen, number of CPUs, memory and OS. This check is 
    useful for parallel processing. 
    
    Returns
    ------------
    system_info: dict
        Dictionary keys refer information about the system 
    """
    print("________________System-Check__________________")
    
    totalMemory = '' # incase system can't figure it out!
    num_threads = ''
    OpSys = ''
    #display processor info
    print("Processor info: %s"%platform.processor())
    num_threads = multiprocessing.cpu_count()
    print("Number of logical CPUs: %i"%num_threads)
    #this message will display if wine is not installed / detected
    helpful_msg ="""   
This version of pyR2 requires wine to run R2.exe, please consider installing
'wine is not an emulator' package @ https://www.winehq.org/. On linux wine can be found on
most reprositories (ubuntu/debian users can use "sudo apt install wine"). Wine acts as
a compatiblity layer between unix like OS systems (ie macOS and linux) and windows programs. 
    """
    msg_flag = False
    #check operating system 
    OpSys=platform.system()    
    if OpSys=='Darwin':
        print("Kernal type: macOS")
    else:
        print("Kernal type: %s"%OpSys)
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
            print("Wine version = "+wine_version)
                          
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
            print("Wine version = "+wine_version)
        except:
            warnings.warn("Wine is not installed!", Warning)
            msg_flag = True
        
    else:
        raise OSError("unrecognised/unsupported operating system")
     
    if totalMemory != '':
        totalMemory = int(totalMemory)
        print("Total RAM available: %i Mb"%totalMemory)
        
        #print some warnings incase the user has a low end PC
        if totalMemory <= 4000:
            warnings.warn("The amount of RAM currently installed is low (<4Gb), complicated ERT problems may incur memory access voilations", Warning)
    
    if num_threads!= '':
        if num_threads <=2:
            warnings.warn("Only one or two CPUs detected, multithreaded workflows will not perform well.", Warning)
            
    if msg_flag:
        print(helpful_msg)
    
    return {'memory':totalMemory,'core_count':num_threads,'OS':OpSys}

#info = systemCheck()
    
