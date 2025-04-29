# GMSH WRAPPER FUNCTIONS
"""
Created on Tue Apr 10 14:32:14 2018 in python 3.6.5
Wrapper for creating 2d and 3d triangular mesh scripts for gmsh and 
parsing meshes into a python environment. 

@author: jimmy Boyd - jamyd91@bgs.ac.uk
"""
#python standard libraries 
import os, warnings
#general 3rd party libraries
import numpy as np
# import matplotlib.pyplot as plt 
import matplotlib.path as mpath
from scipy.spatial import cKDTree

#%% utility functions 
def arange(start,incriment,stop,endpoint=0):#create a list with a range without numpy 
    delta=stop-start
    iterations=int(delta/incriment)
    val=start
    cache=[]
    for i in range(iterations):
        cache.append(val)
        val=val+incriment
    if endpoint==1:
        cache.append(stop)
    return cache

def moving_average(array,N=3):
    """compute the moving average for a list/array. N is the number
    of elements averaged over, it always uses the quiered element as the central
    point and hence must be an odd number""" 
    if not isinstance(N,int) or N%2==0:
        raise Exception("N must be an odd integer")
    idx = round(N/2)-1
    length = len(array)
    out = [0]*length
    for i in range(len(array)):
        if i<idx:#cap minimum index
            a = 0
            b = i + idx + 1
        elif length-i<idx: # cap maximum index
            a = i - idx
            b = len(array)-1
        else:
            a = i - idx
            b = i + idx + 1            
        out[i] = sum(array[a:b])/len(array[a:b])
    return out
    
def find_dist(elec_x, elec_y, elec_z): # find maximum and minimum electrode spacings 
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

#organise unique xy positions clockwise 
def bearing(x,y):
    x0 = 0; y0 = 0 
    if x == x0 and y == y0: # to avoid 0/ error 
        return 0 
    # check for points on bearing and return 
    if x == x0 and y > y0:
        return 0.0
    elif x == x0 and y < y0:
        return 180.0
    elif x > x0 and y==y0:
        return 90.0
    elif x < x0 and y==y0:
        return 270.0 
    # compute raw angle between centre point and point on circle 
    dx = np.abs(x - x0)
    dy = np.abs(y - y0) 
    theta = np.rad2deg(np.arctan(dy/dx)) 
    if x > x0 and y > y0:
        return 90-theta 
    elif x > x0 and y < y0:
        return 90+theta
    elif x < x0 and y < y0:
        return 270-theta  
    elif x < x0 and y > y0:
        return 270+theta 
    
#%% refinement field (depreciated)
def add_ball_field(fh,nfield,x,y,z,r,mincl,maxcl,thick=1):
    fh.write('Field[%i] = Ball;\n'%nfield)

    # set radius and area of influence 
    fh.write('Field[%i].Radius = %f;\n'%(nfield,r))
    fh.write('Field[%i].Thickness = %f;\n'%(nfield,thick))
    # set characteristic lengths 
    fh.write('Field[%i].VIn = %f;\n'%(nfield,mincl))
    fh.write('Field[%i].VOut = %f;\n'%(nfield,maxcl))
    # set ball origin 
    fh.write('Field[%i].XCenter = %f;\n'%(nfield,x))
    fh.write('Field[%i].YCenter = %f;\n'%(nfield,y))
    fh.write('Field[%i].ZCenter = %f;\n'%(nfield,z))
    
def add_box_field(fh, nfield, xmin, xmax, ymin, ymax, zmin, zmax, mincl, maxcl,
                  thick=0):
    # setup box and volume of influence 
    fh.write('Field[%i] = Box;\n'%nfield)
    fh.write('Field[%i].Thickness = %f;\n'%(nfield,thick))
    # set characteristic lengths 
    fh.write('Field[%i].VIn = %f;\n'%(nfield,mincl))
    fh.write('Field[%i].VOut = %f;\n'%(nfield,maxcl))
    # setup geometry 
    fh.write('Field[%i].XMin = %f;\n'%(nfield,xmin))
    fh.write('Field[%i].XMax = %f;\n'%(nfield,xmax))
    fh.write('Field[%i].YMin = %f;\n'%(nfield,ymin))
    fh.write('Field[%i].YMax = %f;\n'%(nfield,ymax))
    fh.write('Field[%i].ZMin = %f;\n'%(nfield,zmin))
    fh.write('Field[%i].ZMax = %f;\n'%(nfield,zmax))
    
    
def set_fields(fh,nfield):
    field_list = [i+1 for i in range(nfield)]
    nfield +=1 
    fh.write('Field[%i] = Max;\n'%nfield)
    field_list_s=str(field_list).replace('[','{').replace(']','}')
    fh.write('Field[%i].FieldsList = %s;\n'%(nfield,field_list_s))
    fh.write('Background Field = %i;\n'%nfield)
    return nfield 

#%% parse a .msh file
def mshParseLegacy(file_path, debug=True):
    if debug: # print outputs? 
        def stream(s,**kwargs):
            print(s,**kwargs)
    else:
        def stream(s):
            pass 
        
    if not isinstance(file_path,str):
        raise Exception("expected a string argument for msh_parser")
    fid=open(file_path,'r')# Open text file
    #Idea: Read Mesh format lines $MeshFormat until $Nodes
    dump=fid.readlines()
    fid.close()
    if len(dump) <= 1:
        raise Exception("Target file is empty!!...aborting!")
    #check the file is a mesh format
    if dump[0].strip() != '$MeshFormat':#removes formating strings, checks if the file is a gmsh file
        raise Exception("Unrecognised file type...aborting!")
    mesh_format=dump[1].strip() 
    formats = ['2.2 0 8','4 0 8']
    if mesh_format not in formats:#warn people that the code was developed with a different mesh format in mind
        warnings.warn('Mesh file format unrecognised ... some errors may occur!\n')
    else:
        stream('Using legacy msh parser...')
    
    stream('Reading %s'%file_path)
    
    if mesh_format == '2.2 0 8':
        stream('Msh file version == 2.x')
        gmshV = 3 # assume its gmsh version 3.06
    else:
        stream('Msh file version == 4.x')
        gmshV = 4 # assume V4 and above 
    
    #find where the nodes start 
    for i, line in enumerate(dump):
        if line.find("$Nodes") != -1:#node flag 
            line = dump[i+1].split()
            no_nodes=int(line[-1])
            node_start = i+2
        elif line.find("$EndNodes") != -1:
            node_end = i
            break # stop the loop, should find the nodes start before the end 
    stream('reading node coordinates...')
    #read in number of nodes - at line 5
    #allocate lists for node numbers and coordinates
    node_num=[0]*no_nodes
    nodex=[0]*no_nodes
    nodey=[0]*no_nodes
    nodez=[0]*no_nodes
    #read in node information
    if gmshV == 3:
        for i in range(node_start, node_end):
            line_info=dump[i].split()
            #convert string info into floats
            line=[float(k) for k in line_info]
            node_idx = int(line[0])
            node_num[node_idx-1]=node_idx
            nodex[node_idx-1]=line[1]
            nodey[node_idx-1]=line[2]
            nodez[node_idx-1]=line[3]
    else:
        c = 0
        i = node_start
        while c <  no_nodes:
            line_info=dump[i].split()
            line = [int(k) for k in line_info] # entity line 
            ent = line[-1] #number of nodes to follow
            ent_start = i #starting line 
            for j in range(ent):
                #print(i+j)
                line_info=dump[ent_start+1+j].split()
                #convert string info into floats
                line = [float(k) for k in line_info]
                node_idx = int(line[0])
                node_num[node_idx-1]=node_idx
                nodex[node_idx-1]=line[1]
                nodey[node_idx-1]=line[2]
                nodez[node_idx-1]=line[3]
                c+=1
                i+=1   
            i+=1
        
    #### read in elements 
    # find where the elements start 
    for i, line in enumerate(dump):
        if line.find("$Elements") != -1:
            element_start = i+2
        if line.find("$EndElements") != -1:
            element_end = i
            break # stop the loop, should find the nodes start before the end 
    
    #engage for loop - this time we want to filter out elements which are not needed 
    nat_elm_num = []#native element number to gmsh
    elm_type = []#element type
    phys_entity = [] # defines the physical entity type the element is assocaited with

    ignored_elements=0#count the number of ignored elements
    
    #determine element type 
    stream('Determining element type...') # this depends a bit on the version of gmsh 
    for i in range(element_start,element_end):
        line = dump[i].split()
        if gmshV == 3:
            elm_type.append(int(line[1]))
        else:
            elm_type.append(len(line)-1)
            
    if gmshV == 3:
        lookin4 = [2,4,6] # number of cell nodes we are looking for 
    else:
        lookin4 = [3,4,6]
    if lookin4[0] in elm_type: # then its triangles 
        npere = 3 # number of nodes per element
    if lookin4[1] in elm_type: # forget triangles its tetrahedra
        npere = 4
    if lookin4[2] in elm_type: # forget triangles its prisms 
        npere = 6  
    
    if npere == 3: 
        stream('Triangle')
        con_matrix = [[],[],[]]
        cell_type = [5]
    elif npere == 4: 
        stream('Tetrahedra')
        con_matrix = [[],[],[],[]]
        cell_type = [10]
    elif npere == 6:
        stream('Prism')
        con_matrix = [[],[],[],[],[],[]]
        cell_type = [13]
    else:
        raise ValueError('Cannot parse mesh becuase the relevant cell types cannot be found')
        
    stream('Reading connection matrix...')
        
    phys = 0
    if gmshV == 3:#parse elements for older gmsh file type 
        for i in range(element_start,element_end):
            splitted = dump[i].split()
            line=[int(k) for k in splitted]
            elmV = line[1]
            if npere == 3:
                elmV+=1 # in this format the flag for triangle elements is 2 so add 1
            #convert string info into ints and cache data
            if npere == elmV:#if line contains number of element vertices we want
                nat_elm_num.append(line[0])
                elm_type.append(line[1]) 
                phys_entity.append(line[4]) 
                for j in range(npere):
                    con_matrix[j].append(line[5+j]-1)
            else:
                ignored_elements += 1
    else: #newer gmsh parser
        c = 0
        i = element_start
        while i < element_end:
            #read in flag line 
            splitted = dump[i].split()
            line=[int(k) for k in splitted]
            nef = line[3] # number of elements to follow 
            elmV = line[2] # number of element vertices 
            phys = line[0] # physical entity 
            if npere == 3:
                elmV+=1 # in this format the flag for triangle elements is 2 so add 1
            if npere == elmV:
                nef_start = i + 1
                for j in range(nef):
                    splitted = dump[nef_start+j].split()
                    line=[int(k) for k in splitted]
                    nat_elm_num.append(line[0])
                    # elm_type.append(line[1]) 
                    phys_entity.append(phys) 
                    for k in range(npere):
                        con_matrix[k].append(line[1+k]-1)
                    c+=1
                    i+=1 
                i+=1
            else:
                ignored_elements += nef
                i += nef+1
            
    stream("ignoring %i elements in the mesh file, as they are not required for R2/R3t"%ignored_elements)
    
    real_no_elements=len(nat_elm_num) #'real' number of elements that we actaully want
    if real_no_elements == 0:
        stream("no elements found... aborting" )
        raise Exception ("No elements have been imported, please check formatting of .msh file")
            
    elm_id = [i+1 for i in range(real_no_elements)]
            
    mesh_dict = {'num_elms':real_no_elements,
                'num_nodes':no_nodes,
                'dump':dump,      
                'node_x':nodex,#x coordinates of nodes 
                'node_y':nodey,#y coordinates of nodes
                'node_z':nodez,#z coordinates of nodes 
                'node_id':node_num,#node id number 
                'elm_id':elm_id,#element id number 
                'node_data':con_matrix,#nodes of element vertices
                'cell_type':cell_type,
                'parameters':phys_entity,#the values of the attributes given to each cell 
                'parameter_title':'regions',
                'dict_type':'mesh_info',
                'original_file_path':file_path} 
    
    stream('Finished reading .msh file')
    
    return mesh_dict # return a python dictionary 


def mshParse47(fname, debug=True):
    if debug: # print outputs? 
        def stream(s,**kwargs):
            print(s,**kwargs)
    else:
        def stream(s,**kwargs):
            pass 
        
    if not isinstance(fname,str):
        raise Exception("expected a string argument for fname")
    fid = open(fname,'r')# Open text file
    #Idea: Read Mesh format lines $MeshFormat until $Nodes
    dump=fid.readlines()
    fid.close()
    if len(dump) <= 1:
        raise Exception("Target file is empty!!...aborting!")
    #check the file is a mesh format
    if dump[0].strip() != '$MeshFormat':#removes formating strings, checks if the file is a gmsh file
        raise Exception("Unrecognised file type...aborting!")
    mesh_format=dump[1].strip() 
    formats = ['4.1 0 8']
    if mesh_format not in formats:#warn people that the code was developed with a different mesh format in mind
        warnings.warn('Mesh file format unrecognised ... some errors may occur!\n') 
    else:
        stream('Reading Msh file version == 4.1')
    
    stream('Reading %s'%fname)
    
    #find where the nodes start
    node_start = -1
    node_end = -1
    for i, line in enumerate(dump):
        if line.find("$Nodes") != -1:#node flag 
            node_start = i
        elif line.find("$EndNodes") != -1:
            node_end = i
            break # stop the loop, should find the nodes start before the end 
    if node_start == -1 or node_end == -1:
        raise ValueError("Could not find $Nodes or $EndNodes in generated .msh file")

    stream('reading node coordinates...')
    #read in node stats 
    line = dump[node_start+1].split()
    node_ent = int(line[0])
    numnp = int(line[3])
    min_tag = int(line[2])
    max_tag = int(line[3])
    
    #allocate arrays for nodes 
    node = np.zeros((numnp,3))
    node_id=[0]*numnp
    #read in node information
    tmp = node_start + 2
    for i in range(node_ent):
        lno = tmp # line number of node block info 
        line = dump[lno].split()
        nodeinblock = int(line[-1])
        for j in range(nodeinblock):
            line_np=dump[lno+1+j].split()#node id line 
            line_co=dump[lno+1+j+nodeinblock].split()#node coordinate line 
            #convert info 
            idx = int(line_np[0])
            node_id[idx-1] = idx
            
            node[idx-1,0]=float(line_co[0])
            node[idx-1,1]=float(line_co[1])
            node[idx-1,2]=float(line_co[2])
        tmp = lno+nodeinblock+nodeinblock+1
        
    #### read in elements 
    # find where the elements start 
    for i, line in enumerate(dump):
        if line.find("$Elements") != -1:
            element_start = i+2
        if line.find("$EndElements") != -1:
            element_end = i
            break # stop the loop, should find the nodes start before the end 
    
    #engage for loop - this time we want to filter out elements which are not needed 
    nat_elm_num = []#native element number to gmsh
    elm_type = []#element type
    phys_entity = [] # defines the physical entity type the element is assocaited with

    ignored_elements=0#count the number of ignored elements
    
    #determine element type 
    stream('Determining element type...',end='') # this depends a bit on the version of gmsh 
    for i in range(element_start,element_end):
        line = dump[i].split()

        elm_type.append(len(line)-1)
            

    lookin4 = [3,4,6]
    if lookin4[0] in elm_type: # then its triangles 
        npere = 3 # number of nodes per element
    if lookin4[1] in elm_type: # forget triangles its tetrahedra
        npere = 4
    if lookin4[2] in elm_type: # forget triangles its prisms 
        npere = 6  
    
    if npere == 3: 
        stream('Triangle')
        con_matrix = [[],[],[]]
        vtk_type = 5
    elif npere == 4: 
        stream('Tetrahedra')
        con_matrix = [[],[],[],[]]
        vtk_type = 10
    elif npere == 6:
        stream('Prism')
        con_matrix = [[],[],[],[],[],[]]
        vtk_type = 13
    else:
        raise ValueError('Cannot parse mesh becuase the relevant cell types cannot be found')
        
    stream('Reading connection matrix...')
        
    phys = 0
    c = 0
    i = element_start
    while i < element_end:
        #read in flag line 
        splitted = dump[i].split()
        line=[int(k) for k in splitted]
        nef = line[3] # number of elements to follow 
        elmV = line[2] # number of element vertices 
        phys = line[1] # physical entity 
        if npere == 3:
            elmV+=1 # in this format the flag for triangle elements is 2 so add 1
        if npere == elmV:
            nef_start = i + 1
            for j in range(nef):
                splitted = dump[nef_start+j].split()
                line=[int(k) for k in splitted]
                nat_elm_num.append(line[0])
                # elm_type.append(line[1]) 
                phys_entity.append(phys) 
                for k in range(npere):
                    con_matrix[k].append(line[1+k]-1)
                c+=1
                i+=1 
            i+=1
        else:
            ignored_elements += nef
            i += nef+1
    # connection = np.array(con_matrix).T
            
    stream("ignoring %i elements in the mesh file, as they are not required for (c)R2/3t"%ignored_elements)
    
    real_no_elements=len(nat_elm_num) #'real' number of elements that we actaully want
    if real_no_elements == 0:
        stream("no elements found... aborting" )
        raise Exception ("No elements have been imported, please check formatting of .msh file")
            
    elm_id = [i+1 for i in range(real_no_elements)]
    cell_type = [vtk_type]*real_no_elements
    
    mesh_dict = {'num_elms':real_no_elements,
                'num_nodes':numnp,
                'dump':dump,      
                'node_x':node[:,0],#x coordinates of nodes 
                'node_y':node[:,1],#y coordinates of nodes
                'node_z':node[:,2],#z coordinates of nodes 
                'node_id':node_id,#node id number 
                'elm_id':elm_id,#element id number 
                'node_data':con_matrix,#nodes of element vertices
                'cell_type':cell_type,
                'parameters':phys_entity,#the values of the attributes given to each cell 
                'parameter_title':'regions',
                'dict_type':'mesh_info',
                'original_file_path':fname} 
    
    stream('Finished reading .msh file')
    
    return mesh_dict # return a python dictionary 

def mshParse(fname, debug=True):
    """Import a gmsh mesh file into ResIPy. 
    
    Parameters
    ----------
    fname: string
        file path to mesh file. note that a error will occur if the file format is not as expected
   
    Returns
    ----------
    Mesh class
    """
    fh = open(fname)
    l0 = fh.readline()
    l1 = fh.readline().strip()
    fh.close()
    if l1 == '4.1 0 8': 
        return mshParse47(fname, debug)
    else:
        return mshParseLegacy(fname, debug)
    
#%% refinment fields 
def addBallField(fh,nfield,x,y,z,r,thick=0,VIn='cl_factor',
                 VOut='cln_factor'):
    fh.write('Field[%i] = Ball;\n'%nfield)

    # set radius and area of influence 
    fh.write('Field[%i].Radius = %f;\n'%(nfield,r))
    fh.write('Field[%i].Thickness = %f;\n'%(nfield,thick))
    # set characteristic lengths 
    fh.write('Field[%i].VIn = cl*%s;\n'%(nfield,VIn))
    fh.write('Field[%i].VOut = cl*%s;\n'%(nfield,VOut))
    # set ball origin 
    fh.write('Field[%i].XCenter = %f;\n'%(nfield,x))
    fh.write('Field[%i].YCenter = %f;\n'%(nfield,y))
    fh.write('Field[%i].ZCenter = %f;\n'%(nfield,z))
    
def setFields(fh,nfield):
    field_list = [i+1 for i in range(nfield)]
    nfield +=1 
    fh.write('Field[%i] = Max;\n'%nfield)
    field_list_s=str(field_list).replace('[','{').replace(']','}')
    fh.write('Field[%i].FieldsList = %s;\n'%(nfield,field_list_s))
    fh.write('Background Field = %i;\n'%nfield)
    return nfield 
        
#%% whole space problems 
def wholespace2d(electrodes, padding = 20, electrode_type = None, geom_input = None,
                 file_path='mesh.geo',cl=-1,cl_factor=50,fmd=None,dp_len=None):
    """Writes a gmsh .geo for a 2D whole space. Ignores the type of electrode. 
    
    Parameters
    ----------
    electrodes: array like
        first column/list is the x coordinates of electrode positions, second column
        is the elevation
    padding: float, optional
        Padding in percent on the size the fine mesh region extent. Must be bigger than 0.
    electrode_type: list
        Electrode types. 
    geom_input: dict, optional
        Allows for further customisation of the 2D mesh, its a
        dictionary contianing surface topography, polygons and boundaries 
    file_path: string, optional 
        name of the generated gmsh file (can include file path also) (optional)
    cl: float, optional
        characteristic length (optional), essentially describes how big the nodes 
        assocaited elements will be. Usually no bigger than 5.
    cl_factor: float, optional
        Describes the growth factor applied to the outer fine mesh region.  
    fmd: N/A
        Parameter is not used. Needed for keyword argument compatiblity with 
        halfspace2d.
    dp_len: N/A
        Parameter is not used. Needed for keyword argument compatiblity with 
        halfspace2d.            
    
    Returns
    ----------
    Node_pos: numpy array
        The indexes for the mesh nodes corresponding to the electrodes input, the ordering of the nodes
        should be the same as the input of 'electrodes'
    .geo: file
        Can be run in gmsh

    NOTES
    ----------
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
    original_x = np.array(electrodes[0])
    original_z = np.array(electrodes[1])   
    elec_x = electrodes[0]
    elec_z = electrodes[1]
    
    if len(elec_x) != len(elec_z):
        raise ValueError("The length of the x coordinate array does not match of the Z coordinate")
    
    if geom_input != None: 
        if not isinstance(geom_input,dict):
            raise TypeError ("'geom_input' is not a dictionary type object. Dict type is expected for the first argument of genGeoFile_adv")
    elif geom_input is None:
        geom_input = {}
        
    if file_path.find('.geo')==-1:
        file_path=file_path+'.geo'#add file extension if not specified already
    
    if cl==-1:
        dist_sort = np.unique(find_dist(elec_x,[0]*len(elec_x),elec_z))
        cl = dist_sort[0]/2 # characteristic length is 1/2 the minimum electrode distance
        
    rem_idx = []
    elec_x_cache = []
    elec_z_cache = []   
    if electrode_type is not None: 
        for i in range(len(electrode_type)):
            if electrode_type[i] == 'remote':
                rem_idx.append(i)
            else:
                elec_x_cache.append(elec_x[i])
                elec_z_cache.append(elec_z[i])
        elec_x_cache = np.array(elec_x_cache)
        elec_z_cache = np.array(elec_z_cache)
    else:
        elec_x_cache = np.array(electrodes[0])
        elec_z_cache = np.array(electrodes[1])
    
    if len(rem_idx) > 0:
        elec_x = np.delete(elec_x,rem_idx)
        elec_z = np.delete(elec_z,rem_idx)
        
    fh = open(file_path,'w') #file handle
    
    fh.write("//Gmsh wrapper code version 1.0 (run the following in gmsh to generate a triangular mesh for 2D whole space)\n")
    fh.write("//2D mesh coordinates\n")
    fh.write("Mesh.Binary = 0;//specify we want ASCII format\n")
    fh.write("cl=%.2f;//define characteristic length\n" %cl)
    
    #create square around all of the electrodes
    x_dist = abs(np.max(elec_x) - np.min(elec_x))
    z_dist = abs(np.max(elec_z) - np.min(elec_z))
    if x_dist<0.2:x_dist=5 # protection against small (or zero) padding 
    if z_dist<0.2:z_dist=5
    max_x = np.max(elec_x) + (padding/100)*x_dist
    min_x = np.min(elec_x) - (padding/100)*x_dist
    max_z = np.max(elec_z) + (padding/100)*z_dist
    min_z = np.min(elec_z) - (padding/100)*z_dist
    
    fh.write("//Fine mesh region.\n")
    fh.write("cl_pad=%f;//padding characteristic length\n"%(cl*cl_factor))#padding applied to fine mesh outer nodes
    #add points to file 
    no_pts = 1
    loop_pt_idx=[no_pts]
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl_pad};\n"%(no_pts, max_x, 0, max_z))
    no_pts += 1
    loop_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl_pad};\n"%(no_pts, max_x, 0, min_z))
    no_pts += 1
    loop_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl_pad};\n"%(no_pts, min_x, 0, min_z))
    no_pts += 1
    loop_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl_pad};\n"%(no_pts, min_x, 0, max_z))
    
    #add line loop
    no_lns = 0 
    for i in range(4):
        no_lns += 1 
        if i == 3:
            fh.write("Line(%i) = {%i,%i};\n"%(no_lns,loop_pt_idx[i],loop_pt_idx[0]))
        else:
            fh.write("Line(%i) = {%i,%i};\n"%(no_lns,loop_pt_idx[i],loop_pt_idx[i+1]))
         
    #Nueman boundary 
    flank_x = 80*x_dist
    flank_z = 50*z_dist 
    fh.write("//Nueman boundary \n")
    cl2 = cl*150
    fh.write("cl2 = %.2f;\n"%cl2)
    no_pts += 1
    nmn_pt_idx=[no_pts]
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl2};\n"%(no_pts, max_x+flank_x, 0, max_z+flank_z))
    no_pts += 1
    nmn_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl2};\n"%(no_pts, max_x+flank_x, 0, min_z-flank_z))
    no_pts += 1
    nmn_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl2};\n"%(no_pts, min_x-flank_x, 0, min_z-flank_z))
    no_pts += 1
    nmn_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl2};\n"%(no_pts, min_x-flank_x, 0, max_z+flank_z))
    
    for i in range(4):
        no_lns += 1 
        if i == 3:
            fh.write("Line(%i) = {%i,%i};\n"%(no_lns,nmn_pt_idx[i],nmn_pt_idx[0]))
        else:
            fh.write("Line(%i) = {%i,%i};\n"%(no_lns,nmn_pt_idx[i],nmn_pt_idx[i+1]))
            
    fh.write("Line Loop(2) = {5,6,7,8};\n")    
    fh.write("Plane Surface(1) = {2};\n")
    fh.write("Line{1,2,3,4} In Surface{1};\n")
    
    fh.write("//Electrode positions.\n")
    node_pos = np.zeros(len(elec_x),dtype=int)
    for i in range(len(elec_x)):
        no_pts += 1
        node_pos[i] = no_pts
        fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl};\n"%(no_pts, elec_x[i], 0, elec_z[i]))
        fh.write("Point{%i} In Surface{1};\n"%(no_pts))# put the point surface
    fh.write("//End electrodes\n")
    
    if len(rem_idx)>0: #then we need to add remote electrodes to the mesh
        fh.write("\n//Remote electrodes \n")
        e_pt_idx = [0]*len(rem_idx)
        #put points in lower left hand corner of the outer region
        #if an electrode does have the remote flag, then chances are it has 
        #a bogus coordinate associated with it (ie. -99999)
        for k in range(len(rem_idx)):
            no_pts += 1
            remote_x = min_x-flank_x + (cl*cl_factor*(k+1))
            remote_z = min_z-flank_z + (cl*cl_factor*(k+1))
            fh.write("Point(%i) = {%.2f,%.2f,%.2f,cln};//remote electrode\n"%(no_pts,remote_x,0,remote_z))
            e_pt_idx[k] = no_pts
            original_x[rem_idx[k]] = remote_x 
            original_z[rem_idx[k]] = remote_z 
            elec_x_cache = np.append(elec_x_cache,original_x[rem_idx[k]])
            elec_z_cache = np.append(elec_z_cache,original_z[rem_idx[k]])
            
        node_pos = np.append(node_pos,e_pt_idx) #add remote electrode nodes to electrode node positions 
        fh.write("Point{%s} In Surface{1};\n"%(str(e_pt_idx).strip('[').strip(']')))
        fh.write('//End of remote electrodes.\n')
    
    fh.write("\n//Adding polygons?\n")
    no_lin=no_lns
    no_plane=0
    count = 0    
    while True:  
        count += 1
        key = 'polygon'+str(count)

        if key in geom_input.keys():
            plyx = geom_input[key][0]
            plyz = geom_input[key][1]
            plyy = [0]*len(plyx)
            pt_idx = [0] *len(plyx)
            fh.write("//polygon vertices for polygon %i\n"%count)
            for k in range(len(plyx)):
                no_pts += 1
                fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl};//polygon vertex \n"%(no_pts,plyx[k],plyy[k],plyz[k]))
                pt_idx[k] = no_pts
            fh.write("//put lines between each vertex\n")
            line_idx = []
            for i in range(len(pt_idx)):
                idx = pt_idx[i]
                no_lin += 1
                if i == len(pt_idx)-1:
                    fh.write("Line (%i) = {%i,%i};\n"%(no_lin,idx,pt_idx[0]))
                else:
                    fh.write("Line (%i) = {%i,%i};\n"%(no_lin,idx,idx+1))
                line_idx.append(no_lin)
            #make line loop out of polygon
            fh.write("//make lines forming polygon into a line loop? - current inactive due to unexpected behaviour in gmsh\n")
            no_lin += 1
            fh.write("//Line Loop(%i) = {%s};\n"%(no_lin,str(line_idx).strip('[').strip(']')))
            no_plane +=1
            fh.write("//Plane Surface(%i) = {%i};\n"%(no_plane,no_lin))
            fh.write("Line{%s} In Surface{1};\n"%str(line_idx).strip('[').strip(']'))
            
        else: 
            fh.write("//end of polygons.\n")
            print('%i polygons added to input file'%(count-1))
            break  

    fh.write("\n//Adding boundaries?\n")
    count = 0   
    while True:
        count += 1        
        key = 'boundary'+str(count)

        if key in geom_input.keys(): 
            bdx = geom_input[key][0]
            bdz = geom_input[key][1]
            bdy = [0]*len(bdx)
            pt_idx = [0] *len(bdx)
            fh.write("// vertices for boundary line %i\n"%count)
            for k in range(len(bdx)):
                no_pts += 1 
                fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl};//boundary vertex \n"%(no_pts,bdx[k],bdy[k],bdz[k]))
                pt_idx[k] = no_pts
            fh.write("//put lines between each vertex\n")
            line_idx = []
            for i in range(len(pt_idx)-1):
                idx = pt_idx[i]
                no_lin += 1
                fh.write("Line (%i) = {%i,%i};\n"%(no_lin,idx,idx+1))
                line_idx.append(no_lin)
            fh.write("Line{%s} In Surface{1};\n"%str(line_idx).strip('[').strip(']'))
                
        else:
            fh.write("//end of boundaries.\n")
            print('%i boundary(ies) added to input file'%(count-1))
            break       
        
    if 'refine' in geom_input.keys():
        fh.write('\n//Refinement points \n')
        px = np.array(geom_input['refine'][0])
        pz = np.array(geom_input['refine'][1])
        if len(geom_input['refine']) == 3: 
            pl = np.array(geom_input['refine'][2])
        else:
            pl = np.ones(px.shape[0])*cl
        #work out which nodes are inside the fine mesh region 
        fmrx = [max_x,max_x,min_x,min_x]
        fmrz = [max_z,min_z,min_z,max_z]
        path = mpath.Path(np.c_[fmrx,fmrz])
        inside = path.contains_points(np.c_[px,pz])
        r_pt_idx = []
        for k in range(len(px)):
            if not inside[k]:
                continue 
            no_pts += 1 
            fh.write("Point(%i) = {%.2f,%.2f,%.2f,%.2f};//refinement point %i\n"%(no_pts, px[k], 0, pz[k], pl[k], k+1))
            r_pt_idx.append(no_pts) 
        fh.write("Point{%s} In Surface{%i};//include nodes in meshing\n"%(str(r_pt_idx).strip('[').strip(']'), 1))
        print('%i refinement points added to input file'%(len(r_pt_idx)))
        fh.write('//End refinement points\n')
    
    fh.close()
    print("writing .geo to file completed, save location:\n%s\n"%os.getcwd())
    
    #sort node ordering back into the original input ordering    
    #find the original indexes  
    original_idx = [0]*len(node_pos)
    for i in range(len(node_pos)):
        idx = (np.abs(elec_x_cache - original_x[i])<1e-15) & (np.abs(elec_z_cache - original_z[i])<1e-15)
        idx = idx.tolist()
        original_idx[i] = idx.index(True)

    ordered_node_pos = node_pos[original_idx].astype(int)
    
    return ordered_node_pos

def wholespace3d(elec_x, elec_y, elec_z = None,
                 fmd=-1, file_path='mesh3d.geo',
                 cl=-1, cl_factor=-1, cln_factor=1000, dp_len=-1, 
                 mesh_refinement=None, use_fields=False, dump=None,
                 flank_fac = 12):
    """
    Writes a gmsh .geo for a 3D wholespace (spherical) mesh. Ignores the type 
    of electrode. Z coordinates should be given as depth below the surface! 
    If Z != 0 then its assumed that the electrode is buried. 
    
    Parameters
    ----------
    elec_x: array like
        electrode x coordinates 
    elec_y: array like 
        electrode y coordinates 
    elec_z: array like 
        electrode z coordinates 
    fmd: float, optional 
        Depth of investigation of the survey. Doesn't actually do anything in this
        function but kept to maintain consistency with half space code.  
    file_path: string, optional 
        name of the generated gmsh file (can include file path also) (optional)
    cl: float, optional
        characteristic length (optional) of electrodes, essentially describes how big the nodes 
        assocaited elements will be on the electrodes. Usually no bigger than 5. If set as -1 (default)
        a characteristic length 1/4 the minimum electrode spacing is computed.
    cl_factor: float, optional 
        Characteristic length mulitplier for the sub-surface points on away from the 
        electrodes on the top of the fine mesh region.  This allows for tuning of 
        the incremental size increase with depth in the mesh, usually set to 2 such 
        that the elements at the DOI are twice as big as those
        at the surface. The reasoning for this is because the sensitivity of ERT drops
        off with depth. 
    cln_factor: float, optional
        Characteristic length mulitplier for the Neumann boundary points in the 
        coarse mesh region. 
    mesh_refinement: dict, pd.DataFrame, optional 
        Coordinates for discrete points in the mesh (advanced use cases). 
    dump : function, optional
        If None, output is printed using `print()`. Else the given function is passed.
    flank_fac: float, int, optional
        Defines how far away the outer radius of the of the mesh. The radius is 
        flank_factor * dp_len. 
    
    Returns
    ----------
    Node_pos: numpy array
        The indexes for the mesh nodes corresponding to the electrodes input, the ordering of the nodes
        should be the same as the input of the electrodes. 
    """
    
    if dump is None:
        def dump(x):
            print(x)
    if dp_len != -1 and dp_len<0:
        raise ValueError('Can not have a negative dipole length')
    if fmd != -1 and fmd<0:#then set to a default 
        raise ValueError('Can not have a negative depth of fine mesh')
    if cl != -1 and cl<0:
        raise ValueError('Can not have a negative mesh characteristic length')
        
    if file_path.find('.geo')==-1:
        file_path=file_path+'.geo'#add file extension if not specified already
    
    # check for z coordinates 
    if elec_z is None:
        elec_z = [0]*len(elec_x)

    # definition variable setup  
    nelec = len(elec_x)
    no_pts = 0
    nfield = 0 
    node_pos = np.zeros(nelec,dtype=int)
    dist = np.unique(find_dist(elec_x, elec_y, elec_z)) 
    
    dmin = dist[1]
    dmax = dist[-1]
    xmean = np.mean(elec_x)
    ymean = np.mean(elec_y)
    
    if cl==-1:
        cl = dmin/4 # characteristic length is 1/4 the minimum electrode distance

    if cl_factor == -1: 
        cl_factor = 5 
        
    if dp_len == -1: # compute largest possible dipole length 
        dp_len = dmax # maximum possible dipole length
    
    fmd = None 
    outer_rad = dp_len*flank_fac # outer raduis of hemi-sphere 
    
    # create file 
    fh = open(file_path,'w')
    fh.write('//3D Wholespace (tetra) mesh for ResIPy\n')
    fh.write("Mesh.Binary = 0;//specify we want ASCII format\n")
    fh.write('SetFactory("OpenCASCADE");\n')
    fh.write('cl=%f;//electrode characteristic length\n'%float(cl))
    fh.write('cl_factor=%f;//Fine region characteristic factor\n'%float(cl_factor))
    fh.write('cln_factor=%f;//Neumann characteristic factor\n'%float(cln_factor))
    
    template = 'Sphere(1) = {%f, %f, %f, %f, -Pi/2, Pi/2, Pi*2};//create a sphere\n' 
    no_pts+=2 # seems adding a sphere adds in 2 points 
    fh.write(template%(xmean,ymean,0,outer_rad))#
    fh.write("Mesh.CharacteristicLengthFromPoints=1;\n")
    # fh.write("MeshSize{1,2} = cl*cln_factor;")
    
    # create electrodes
    fh.write('//Start electrodes, include refinement sphere around the electrode\n')
    template = "Point(%i) = {%f, %f, %f, cl};//electrode coordinate\n"
    for i in range(nelec):
        no_pts += 1 
        fh.write(template%(no_pts,elec_x[i],elec_y[i],elec_z[i]))
        node_pos[i] = no_pts
        fh.write("Point{%i} In Volume{1};//specify electrode volume\n"%no_pts)
        if use_fields: 
            nfield += 1 
            addBallField(fh,nfield,elec_x[i],elec_y[i],elec_z[i],
                           dmin/4,0)
    if use_fields:
        set_fields(fh, nfield)
    fh.write("//End of electrodes\n")
        
    # check if any mesh refinement is requested 
    if mesh_refinement is not None:
        fh.write('//Custom refinement points\n')
        # find surface points 
        rx = mesh_refinement['x']
        ry = mesh_refinement['y']
        rz = mesh_refinement['z']
        npoints = len(rx) # number of refinement points 
        if 'cl' in mesh_refinement.keys():
            rcl = mesh_refinement['cl'] # use cl specified for each point 
            fh.write('//Refinement point characteristic length specified for each point individually\n')
        else:
            rcl = [cl]*npoints # use same cl as for electrodes 
            fh.write('//Refinement point characteristic length the same as for electrodes\n')
        for i in range(npoints):
            no_pts += 1
            fh.write("Point (%i) = {%.16f, %.16f, %.16f, %.16f};\n"%(no_pts, rx[i], ry[i], rz[i], rcl[i]))
            fh.write("Point{%i} In Volume{%i};//buried refinement point\n"%(no_pts,1))# put the point in volume 
        fh.write('//End mesh refinement points\n')
        
    fh.write("//End of gmsh script\n")
        
    fh.close()
    
    return node_pos 


#%% 2D halfspace 
def halfspace2d(electrodes, electrode_type = None, geom_input = None,
                file_path='mesh.geo',fmd=-1,dp_len=-1,cl=-1,cl_factor=2.0,
                edge_factor=5, debug=False):
    """Writes a gmsh .geo file for a 2d study area with topography assuming we wish to add electrode positions
    
    Parameters
    ----------
    electrodes : array like
        2 Column/list. First column is the x coordinates of electrode positions, second
        is the elevation of electrode positions.
    electrode_type : list, optional
        List should be the same length as the electrode coordinate argument. Each entry in
        the list references the type of electrode: 
            'electrode' = surface electrode coordinate, will be used to construct the topography in the mesh
            'buried' = buried electrode, placed the mesh surface
            'borehole' = borehole electrode, electrodes will be placed in the mesh with a line connecting them. 
                        borehole numbering starts at 1 and ascends numerically by 1. 
            'remote' = remote electrode, this will not be actaully meshed, as often remote electrode coordinates
                        are given arbitary values like -999. Instead every remote electrode found will be 
                        randomly added as a node in the background region towards the lower left hand corner of 
                        of the mesh. 
    geom_input : dict, optional
        Allows for further customisation of the 2D mesh, its a
        dictionary contianing surface topography, polygons and boundaries 
    file_path : string, optional 
        Name of the generated gmsh file (can include file path also) (optional)
    fmd : float, optional 
        defines depth of fine mesh region, this should be depth below ground 
    cl : float, optional
        Characteristic length (optional), essentially describes how big the nodes 
        assocaited elements will be. Usually no bigger than 5. 
    cl_factor : float, optional 
        Characteristic length factor, this allows for tuning of the incrimental 
        size increase with depth in the mesh, usually set to 2 such that the 
        elements at the DOI are twice as big as those at the surface. The reasoning
        for this is because the sensitivity of ERT drops off with depth. 
    dp_len : float, optional 
        Largest dipole length in the 2D array. Default is the largest electrode 
        spacing. Controls the multipier applied to the edges of the coarse mesh region. 
    edge_factor: float, optional 
        Edge of the mesh is edge_factor*dp_len from the edge of the fine mesh zone. 
        Normally a value of 5 is sufficient. 
    debug : bool, optional
        If `True`, debug messages will be displayed.
    
    Returns
    ----------
    Node_pos: numpy array
        The indexes for the mesh nodes corresponding to the electrodes input, the ordering of the nodes
        should be the same as the input of 'electrodes'
    .geo: file
        Can be run in gmsh

    NOTES
    ----------
     geom_input format:
        the code will cycle through numerically ordered keys (strings referencing objects in a dictionary"),
        currently the code expects a 'surface' and 'electrode' key for surface points and electrodes.
        the first borehole string should be given the key 'borehole1' and so on. The code stops
        searching for more keys when it cant find the next numeric key. Same concept goes for adding boundaries
        and polygons to the mesh. A boundary is a 1D line which is considered to span the fine mesh region,
        and a polygon is a discrete shape within the fine mesh region. See below example for format:
            
            geom_input = {'surface': [surf_x,surf_z],
              'boundary1':[bound1x,bound1z],
              'polygon1':[poly1x,poly1z],
              'refine':[pointx,pointz]}
            
    electrodes and electrode_type (if not None) format: 
        
            electrodes = [[x1,x2,x3,...],[y1,y2,y3,...]]
            electrode_type = ['electrode','electrode','buried',...]
        
        like with geom_input, boreholes should be labelled borehole1, borehole2 and so on.
        The code will cycle through each borehole and internally sort them and add them to 
        the mesh. Surface electrodes can have the flag 'surface' or 'electrode', independently 
        buried electrodes have the tag 'buried'. 
        
    The code expects that all polygons, boundaries and electrodes fall within x values 
    of the actaul survey area. So make sure your topography / surface electrode points cover 
    the area you are surveying, otherwise some funky errors will occur in the mesh. 
    """
    if dp_len != -1 and dp_len<0:
        raise ValueError('Can not have a negative dipole length')
    if fmd != -1 and fmd<0:#then set to a default 
        raise ValueError('Can not have a negative depth of fine mesh')
    if cl != -1 and cl<0:
        raise ValueError('Can not have a negative mesh characteristic length')
        
    if debug:
        def dump(x):
            print(x)
    else:
        def dump(x):
            pass
        
    dump('Generating gmsh input file...\n')
    #formalities and error checks
    if geom_input is not None: 
        if not isinstance(geom_input,dict):
            raise TypeError ("'geom_input' is not a dictionary type object. Dict type is expected for the first argument of genGeoFile")
    else:
        geom_input = {}
    if len(electrodes[0])!=len(electrodes[1]):
        raise ValueError('The length of the electrode x and z arrays does not match')
        
    original_x = np.array(electrodes[0])
    original_z = np.array(electrodes[1])
    bh_flag = False
    bu_flag = False
    #determine the relevant node ordering for the surface electrodes? 
    if electrode_type is not None:
        if not isinstance(electrode_type,list):
            raise TypeError("electrode_type variable is given but expected a list type argument")
        if len(electrode_type) != len(electrodes[0]):
            raise ValueError("mis-match in length between the electrode type and number of electrodes")
        
        surface_idx=[]#surface electrode index
        bh_idx=[]#borehole index
        bur_idx=[]#buried electrode index 
        rem_idx=[]#remote electrode index
        for i,key in enumerate(electrode_type):
            if key == 'electrode': surface_idx.append(i)
            if key == 'surface': surface_idx.append(i)
            if key == 'buried': bur_idx.append(i); bu_flag = True
            if key == 'borehole1': bh_flag = True
            if key == 'remote': rem_idx.append(i)

        if len(surface_idx)>0:# then surface electrodes are present
            elec_x=np.array(electrodes[0])[surface_idx]
            elec_z=np.array(electrodes[1])[surface_idx]
            
        elif len(surface_idx)==0 and 'surface' not in geom_input:
            warnings.warn('All electrodes are buried but no surface topography has been provided, this problem may be better suited as a whole space')
            #fail safe if no surface electrodes are present to generate surface topography 
            if not isinstance(geom_input, dict):
                geom_input = {}
            max_idx = np.argmax(electrodes[0])
            min_idx = np.argmin(electrodes[0])
            topo_x = [electrodes[0][min_idx],electrodes[0][max_idx]]
            z_idx = np.argmax(electrodes[1])
            topo_z = [electrodes[1][z_idx]+1,electrodes[1][z_idx]+1] # add one so we don't cut into the buried in electrodes with the mesh
            geom_input['surface'] = [topo_x,topo_z]
            elec_x = np.array([])
            elec_z = np.array([])

        elif len(surface_idx)==0 and 'surface' in geom_input:
            elec_x = np.array([])
            elec_z = np.array([])
    
    else:
        rem_idx=[]#remote electrode index
        elec_x = np.array(electrodes[0])
        elec_z = np.array(electrodes[1])

    if 'surface' not in geom_input: # make extra topography points if there are none 
        min_idx = np.argmin(elec_x)
        max_idx = np.argmax(elec_x)   
        topo_x = [elec_x[min_idx] - 5*np.mean(np.diff(np.unique(elec_x))),
                  elec_x[max_idx] + 5*np.mean(np.diff(np.unique(elec_x)))]
        topo_z = [elec_z[min_idx],elec_z[max_idx]]
        if bh_flag or bu_flag: 
            distances = find_dist(electrodes[0], np.zeros_like(electrodes[0]), electrodes[1])
            dist_sort = np.unique(distances)
            elecspacing = dist_sort[1]
            topo_x = [elec_x[min_idx] - 2*elecspacing,
                      elec_x[max_idx] + 2*elecspacing]
    else:
        topo_x = geom_input['surface'][0]
        topo_z = geom_input['surface'][1]

    #catch buried electrode special cases 
    if bh_flag or bu_flag: # special case 
        tmp_x = np.array(electrodes[0])
        tmp_y = np.zeros_like(tmp_x)
        tmp_z = np.array(electrodes[1])
        idx = np.ones_like(tmp_x,dtype=bool)
        if len(rem_idx)>0:
            idx[rem_idx] = False
        distances = find_dist(tmp_x[idx], tmp_y[idx], tmp_z[idx])
        dist_sort = np.unique(distances)
        elecspacing = dist_sort[1]
        if dp_len == -1:#there is no change dipole length, maybe due to 1 x coordinate 
            dp_len = dist_sort[-1]
        if fmd == -1:
            fmd = (max(topo_z) - min(tmp_z))+(5*elecspacing)
        if cl == -1:
            cl = elecspacing/4
                  
    if dp_len == -1 and len(elec_x)>0:#compute maximum dipole length
        dp_len = abs(np.max(elec_x) - np.min(elec_x))
    if dp_len == 0: #if its still 0 
        dp_len = 5 # insert obligitory value here
    if fmd == -1:#then set to a default 
        #fmd = np.max(electrodes[1]) - (np.min(electrodes[1]) - abs(np.max(electrodes[0]) - np.min(electrodes[0]))/2)
        fmd = dp_len/3

    dump("fmd in gmshWrap.py: %f"%fmd)
    dump("dp_len in gmshWrap.py: %f"%dp_len)
    
    if cl == -1:
        cl = abs(np.mean(np.diff(elec_x))/4)
            
    if len(topo_x) != len(topo_z):
        raise ValueError("topography x and z arrays are not the same length!")
    if len(elec_x) != len(elec_z):
        raise ValueError("electrode x and z arrays are not the same length!")
        
    if file_path.find('.geo')==-1:
        file_path=file_path+'.geo'#add file extension if not specified already
    
    #start to write the file  
    fh = open(file_path,'w') #file handle
    
    fh.write("//2D mesh script for ResIPy (run the following in gmsh to generate a triangular mesh with topograpghy)\n")
    fh.write("Mesh.Binary = 0;//specify we want ASCII format\n")
    fh.write("cl=%.2f;//define characteristic length\n" %cl)
    fh.write("cl_factor=%2.f;//define characteristic length factor\n"%cl_factor)
    fh.write("//Define surface points\n")
    #we have surface topograpghy, and electrode positions to make use of here:
    x_pts=np.append(topo_x,elec_x)#append our electrode positions to our topograpghy 
    z_pts=np.append(topo_z,elec_z)
    flag=['topography point']*len(topo_x)
    flag=flag+(['electrode']*len(elec_x))   

    # check that all non remote electrodes fall within the fine mesh region 
    elec_x_check = np.array(electrodes[0])
    if len(rem_idx) > 0:
        elec_x_check = np.delete(elec_x_check, rem_idx)
    
    #deal with end case electrodes, check max topo points are outside survey bounds 
    while min(elec_x_check) <= min(x_pts):
        min_idx = np.argmin(x_pts)
        x_pts = np.append(x_pts,x_pts[min_idx] - (cl*cl_factor)) # in this case extend the survey bounds beyond the first (leftmost) electrode 
        z_pts = np.append(z_pts,z_pts[min_idx])
        flag.append('topography point')#add a flag
        
    while max(elec_x_check) >= max(x_pts):
        max_idx = np.argmax(x_pts)
        x_pts = np.append(x_pts,x_pts[max_idx] + (cl*cl_factor)) # in this case extend the survey bounds beyond the last (rightmost) electrode 
        z_pts = np.append(z_pts,z_pts[max_idx])
        flag.append('topography point')
            
    idx=np.argsort(x_pts) # now resort the x positions in ascending order
    x_pts=x_pts[idx] # compute sorted arrays
    z_pts=z_pts[idx]

    #add flags which distinguish what each point is 
    flag_sort=[flag[i] for i in idx]
    
    #we need protection against repeated points, as this will throw up an error in R2 when it comes calculating element areas
    pts = np.c_[x_pts,z_pts]
    tree = cKDTree(pts)
    dist,neigh = tree.query(pts,2)
     
    if any(dist[:,1]<1e-15): # raise an issue when repeated nodes are present 
        # warnings.warn("Duplicated surface and electrode positions were detected, R2 inversion likley to fail due to the presence of elements with 0 area.",Warning)
        #if duplicated points were dectected we should remove them??
        pblm_idx = np.arange(pts.shape[0])[dist[:,1] < 1e-15]
        
        ## go through pairs of repeated nodes and keep one of them, so pair up the repeated nodes 
        pairs = [] # stores pairs of repeated points 
        probed = [] # chaches any points already probed (hence can be ignored)
        del_idx = [] # indexes of points to be deleted 
        for i in pblm_idx:
            pair = (i,neigh[i,1])
            trigger = False # triggers if pairing not already in list 
            if pair[0] not in probed:
                probed.append(pair[0])
                trigger = True 
            if pair[1] not in probed:
                probed.append(pair[1])
                trigger = True 
            if trigger: 
                pairs.append(pair)
                if flag_sort[pair[0]] == flag_sort[pair[1]]: # if flags of same type keep the first 
                    del_idx.append(pair[1])
                else: # keep the electrode point 
                    for p in pair:
                        if flag_sort[p] != 'electrode':
                            del_idx.append(p)

        ## print out to user the problematic nodes 
        error = 'The following surface nodes are repeated!\n'
        # print(error.strip())
        for pair in pairs:
            for p in pair: 
                x_ = x_pts[p]
                z_ = z_pts[p]
                f_ = flag_sort[p]
                line = "X = {:16.8f}, Z = {:16.8f}, flag = {:}\n".format(x_,z_,f_)
                # print(line.strip())
                error += line 
        
        fh.close()
        raise ValueError(error)
        
        #now delete entries from relevant arrays (commented out for now as it cuases issues elsewhere)
        # print('Warning: %i point(s) deleted to avoid meshing artefacts!'%len(del_idx))
        # flag_sort = np.delete(flag_sort,np.array(del_idx))
        # x_pts = np.delete(x_pts ,np.array(del_idx))
        # z_pts = np.delete(z_pts ,np.array(del_idx))
        
    elec_x_cache = x_pts[np.array(flag_sort)=='electrode']
    elec_z_cache = z_pts[np.array(flag_sort)=='electrode']    
        
    
    #now add the surface points to the file
    dump('adding surface points and electrodes to input file...')
    tot_pnts=0#setup a rolling total for points numbering
    sur_pnt_cache=[] # surface points cache 
    for i in range(len(x_pts)):
        tot_pnts+=1
        if i==0 or i==len(x_pts)-1:
            fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl*cl_factor};//%s\n"%(tot_pnts,x_pts[i],0,z_pts[i],flag_sort[i]))
        else:
            fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl};//%s\n"%(tot_pnts,x_pts[i],0,z_pts[i],flag_sort[i]))
        sur_pnt_cache.append(tot_pnts)
    
    #make the lines between each point
    fh.write("//construct lines between each surface point\n")
    tot_lins=0
    sur_ln_cache = []
    for i in range(len(x_pts)-1):
        tot_lins=tot_lins+1
        fh.write("Line(%i) = {%i,%i};\n"%(tot_lins,sur_pnt_cache[i],sur_pnt_cache[i+1]))
        sur_ln_cache.append(tot_lins)
        
    fh.write("//add points below surface to make a fine mesh region\n")#okay so we want to add in the lines which make up the base of the survey area
    if bh_flag: #not allowing mesh refinement for boreholes currently
        cl_factor = 1   
    
    #reflect surface topography at base of fine mesh area, should be a smoothed version of the topography
    extreme_check = np.max(np.diff(z_pts))>3*np.mean(np.diff(x_pts))
    if extreme_check or len(x_pts)<5: # catched few x points or complicated topography 
        div = 1
    else:
        div = 3
    x_base = np.linspace(x_pts[0],x_pts[-1],int(len(x_pts)/div))
    z_base = np.interp(x_base,x_pts,z_pts) - abs(fmd)
    
    basal_pnt_cache = []
    for i in range(len(x_base)):
        tot_pnts += 1
        fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl*cl_factor};//%s\n"%(tot_pnts,x_base[i],0,z_base[i],
                 'base of smoothed mesh region'))
        basal_pnt_cache.append(tot_pnts)
    
    fh.write("//make lines between base of fine mesh region points\n")
    basal_ln_cache=[]
    for i in range(len(x_base)-1):
        tot_lins += 1
        fh.write("Line(%i) = {%i,%i};\n"%(tot_lins,basal_pnt_cache[i],basal_pnt_cache[i+1]))
        basal_ln_cache.append(tot_lins)
    
    fh.write("\n//Adding boundaries\n") # boundaries which stretch across the mesh
    count = 0   
    end_pt_l = [sur_pnt_cache[0]]
    end_pt_r = [sur_pnt_cache[-1]]
    blines = [] # boundary line cache 
    m = (cl_factor-1)/abs(fmd)
    while True:
        count += 1        
        key = 'boundary'+str(count)
        if key in geom_input.keys():
            bdx = np.array(geom_input[key][0])
            bdz = np.array(geom_input[key][1])
            #filtx = (bdx > x_pts[0]) & (bdx < x_pts[-1]) # filter out points beyond x bounds? 
            pt_idx = [0] * (len(bdx)+2)
            # interpolate right and left boundary points 
            x_end = np.array([x_pts[0],x_pts[-1]])
            z_end = np.interp(x_end,bdx,bdz)
            bdx = np.append(bdx,x_end)
            bdz = np.append(bdz,z_end)
            zdepth = np.abs(np.interp(bdx,x_pts,z_pts) - bdz)
            clb = zdepth*m + 1 # boundary characteristic length 
            
            sortid = np.argsort(bdx)
            bdx = bdx[sortid] # reassign to sorted arrays by x coordinate 
            bdz = bdz[sortid] # reassign to sorted arrays by x coordinate 
            
            fh.write("// vertices for boundary line %i\n"%count)
            for k in range(len(bdx)):
                tot_pnts += 1 
                fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl*%.2f};//boundary vertex \n"%(tot_pnts,bdx[k],0,bdz[k],clb[k]))
                pt_idx[k] = tot_pnts
            fh.write("// lines for boundary line %i\n"%count)
            line_idx = []
            end_pt_l.append(pt_idx[0])
            end_pt_r.append(pt_idx[-1])
            
            for i in range(len(pt_idx)-1):
                idx = pt_idx[i]
                tot_lins += 1
                fh.write("Line (%i) = {%i,%i};\n"%(tot_lins,idx,idx+1))
                blines.append(tot_lins)
                
        else:
            fh.write("//end of boundaries.\n")
            dump('%i boundary(ies) added to input file'%(count-1))
            break  
        
    end_pt_l.append(basal_pnt_cache[0])
    end_pt_r.append(basal_pnt_cache[-1])
    end_ln_l = []
    end_ln_r = []
    
    #add left and rightmost lines to close off the fine mesh region
    fh.write("//Add lines at leftmost side of fine mesh region.\n")
    for i in range(len(end_pt_l)-1):
        tot_lins=tot_lins+1;# add to the number of lines rolling total 
        fh.write("Line(%i) = {%i,%i};\n"%(tot_lins,end_pt_l[i],end_pt_l[i+1]))#line from first point on surface to depth
        end_ln_l.append(tot_lins)
    
    fh.write("//Add lines at rightmost side of fine mesh region.\n")
    for i in range(len(end_pt_r)-1):
        tot_lins=tot_lins+1;# add to the number of lines rolling total 
        fh.write("Line(%i) = {%i,%i};\n"%(tot_lins,end_pt_r[i],end_pt_r[i+1]))#line from first point on surface to depth
        end_ln_r.append(tot_lins)

    
    #compile line numbers into a line loop.
    fh.write("//compile lines into a line loop for a mesh surface/region.\n")
    sur_ln_cache_flipped = list(np.flipud(np.array(sur_ln_cache))*-1)
    fine_msh_loop = end_ln_l + basal_ln_cache + [end_ln_r[i]*-1 for i in range(len(end_ln_r))] + sur_ln_cache_flipped
    fh.write("Line Loop(1) = {%s};\n"%str(fine_msh_loop).strip('[').strip(']')) # line loop for fine mesh region 
    
    #now extend boundaries beyond flanks of survey area (so generate your Neummon boundary)
    fh.write("\n//Background region (Neumann boundary) points\n")
    cl_factor2=25*cl_factor#characteristic length multipleier for Neumann boundary 
    cln=cl*cl_factor2#assign new cl, this is so mesh elements get larger from the main model
    fh.write("cln=%.2f;//characteristic length for background region\n" %cln)
    #Background region propeties, follow rule of thumb that background should be 5*largest dipole 
    flank=5*dp_len
    b_max_depth=-abs(fmd)-(3*dp_len)#background max depth
    #add Neumann boundaries on left hand side
    n_pnt_cache=[0,0,0,0]#cache for the indexes of the neumon boundary points 
    tot_pnts=tot_pnts+1;n_pnt_cache[0]=tot_pnts
    fh.write("Point(%i) = {%.2f,%.2f,%.2f,cln};//far left upper point\n"%(tot_pnts,x_pts[0]-flank,0,z_pts[0]))
    tot_pnts=tot_pnts+1;n_pnt_cache[1]=tot_pnts
    fh.write("Point(%i) = {%.2f,%.2f,%.2f,cln};//far left lower point\n"%(tot_pnts,x_pts[0]-flank,0,b_max_depth))
    #add Neumann boundary points on right hand side
    tot_pnts=tot_pnts+1;n_pnt_cache[2]=tot_pnts
    fh.write("Point(%i) = {%.2f,%.2f,%.2f,cln};//far right upper point\n"%(tot_pnts,x_pts[-1]+flank,0,z_pts[-1]))
    tot_pnts=tot_pnts+1;n_pnt_cache[3]=tot_pnts
    fh.write("Point(%i) = {%.2f,%.2f,%.2f,cln};//far right lower point\n"%(tot_pnts,x_pts[-1]+flank,0,b_max_depth))
    #make lines encompassing all the points - counter clock wise fashion
    fh.write("//make lines encompassing all the background points - counter clock wise fashion\n")
    
    n_ln_cache=[0,0,0,0,0]
    tot_lins=tot_lins+1;n_ln_cache[0]=tot_lins
    fh.write("Line(%i) = {%i,%i};\n"%(tot_lins,sur_pnt_cache[0],n_pnt_cache[0]))
    tot_lins=tot_lins+1;n_ln_cache[1]=tot_lins
    fh.write("Line(%i) = {%i,%i};\n"%(tot_lins,n_pnt_cache[0],n_pnt_cache[1]))
    tot_lins=tot_lins+1;n_ln_cache[2]=tot_lins
    fh.write("Line(%i) = {%i,%i};\n"%(tot_lins,n_pnt_cache[1],n_pnt_cache[3]))
    tot_lins=tot_lins+1;n_ln_cache[3]=tot_lins
    fh.write("Line(%i) = {%i,%i};\n"%(tot_lins,n_pnt_cache[3],n_pnt_cache[2]))
    tot_lins=tot_lins+1;n_ln_cache[4]=tot_lins
    fh.write("Line(%i) = {%i,%i};\n"%(tot_lins,n_pnt_cache[2],sur_pnt_cache[-1]))
    
    fh.write("//Add line loops and plane surfaces for the Neumann region\n")
    #now add background region line loop (cos this be made more efficent?)
    basal_ln_cache_flipped = list(np.flipud(np.array(basal_ln_cache))*-1)
    coarse_msh_loop = n_ln_cache + end_ln_r + basal_ln_cache_flipped + [int(end_ln_l[i])*-1 for i in range(len(end_ln_l))]
    fh.write("Line Loop(2) = {%s};\n"%str(coarse_msh_loop).strip('[').strip(']')) # line loop for fine mesh region 
    fh.write("Plane Surface(1) = {1, 2};//Coarse mesh region surface\n")
    
    #now we want to return the point values of the electrodes, as gmsh will assign node numbers to points
    #already specified in the .geo file. This will needed for specifying electrode locations in R2.in   
    node_pos=[i+1 for i, j in enumerate(flag_sort) if j == 'electrode']
    node_pos = np.array(node_pos)
    
    #add borehole vertices and line segments to the survey mesh
    
    no_lin=tot_lins#+1
    no_pts=tot_pnts#+1
    #add borehole electrode strings
    #dump('probing for boundaries and other additions to the mesh')
    count = 0
    if bh_flag:
        fh.write("\n//Boreholes \n")
        while True:
            count += 1
            key = 'borehole'+str(count)#count through the borehole keys
            bh_idx=[]#borehole index
            for i,entry in enumerate(electrode_type):
                if entry == key: bh_idx.append(i)
            #break if no entries found
            if len(bh_idx)==0:
                fh.write("//End of borehole electrodes\n")
                break
                
            bhx = np.array(electrodes[0])[bh_idx]#get borehole coordinate information
            bhz = np.array(electrodes[1])[bh_idx]
            #cache the x and y coordinates 
            elec_x_cache = np.append(elec_x_cache,bhx)
            elec_z_cache = np.append(elec_z_cache,bhz)
            e_pt_idx = [0] *len(bhx)
            fh.write("// string electrodes for borehole %i\n"%count)
            for k in range(len(bhx)):
                no_pts += 1 
                fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl};//borehole %i electrode\n"%(no_pts,bhx[k],0,bhz[k],count))
                e_pt_idx[k] = no_pts
            
            fh.write("//put lines between each electrode\n")
            line_idx = []
            node_pos = np.append(node_pos,e_pt_idx) #add borehole nodes to electrode node positions 
            for i in range(len(e_pt_idx)-1):
                idx = e_pt_idx[i]
                no_lin += 1
                fh.write("Line (%i) = {%i,%i};//borehole %i segment\n"%(no_lin,idx,idx+1,count))
                line_idx.append(no_lin)
            
            fh.write("Line{%s} In Surface{2};\n"%str(line_idx).strip('[').strip(']'))
    
    no_plane = 1 # number of plane surfaces so far (actually two)
    fh.write("\n//Adding polygons\n")
    line_loops = []
    count = 0    
    while True:  
        count += 1
        key = 'polygon'+str(count)

        if key in geom_input.keys():
            plyx = geom_input[key][0]
            plyz = geom_input[key][1]
            plyy= [0]*len(plyx)
            
            pt_idx = [0] *len(plyx)
            fh.write("//polygon vertices for polygon %i\n"%count)
            for k in range(len(plyx)):
                no_pts += 1
                fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl};//polygon vertex \n"%(no_pts,plyx[k],plyy[k],plyz[k]))
                pt_idx[k] = no_pts
            fh.write("//put lines between each vertex\n")
            line_idx = []
            for i in range(len(pt_idx)):
                idx = pt_idx[i]
                no_lin += 1
                if i == len(pt_idx)-1:
                    fh.write("Line (%i) = {%i,%i};\n"%(no_lin,idx,pt_idx[0]))
                else:
                    fh.write("Line (%i) = {%i,%i};\n"%(no_lin,idx,idx+1))
                line_idx.append(no_lin)
            #make line loop out of polygon
            fh.write("//make lines forming polygons from line loops\n")
            no_lin += 1
            fh.write("Line Loop(%i) = {%s};\n"%(no_lin,str(line_idx).strip('[').strip(']')))
            line_loops.append(no_lin)
            no_plane +=1
            fh.write("Plane Surface(%i) = {%i};\n"%(no_plane,no_lin))
            fh.write("Physical Surface (%i) = {%i};\n"%(no_plane,no_plane))
            
        else: 
            fh.write("//end of polygons.\n")
            dump('%i polygons added to input file'%(count-1))
            break  
        
    line_loops = np.array(line_loops)
    no_plane += 1
    if len(line_loops) > 0:
        fh.write("Plane Surface(%i) = {%s, 1};//Fine mesh region surface\n"%(no_plane, ', '.join(line_loops.astype(str))))
    else:
        fh.write("Plane Surface(%i) = {1};//Fine mesh region surface\n"%(no_plane))
    fh.write("\n//Make a physical surface\n")
    fh.write("Physical Surface(1) = {%i, 1};\n"%no_plane)
    
    if len(blines)>0:#add boundary lines if present 
        fh.write("Line{%s} In Surface{2};//boundary lines\n"%str(blines).strip('[').strip(']'))
    
    # add buried electrodes? (added after as we need Surface 1 to be defined)      
    if bu_flag:
        nfield = 0 
        fh.write("\n//Buried electrodes \n")  
        buried_x = np.array(electrodes[0])[bur_idx]#get buried electrode coordinate information
        buried_z = np.array(electrodes[1])[bur_idx]
        buried_y = [0]*len(buried_x)
        e_pt_idx = [0]*len(buried_x)
        for k in range(len(buried_x)):
            no_pts += 1 
            nfield += 1 
            fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl};//buried electrode %i\n"%(no_pts,buried_x[k],buried_y[k],buried_z[k],k+1))
            e_pt_idx[k] = no_pts
            # add_ball_field(fh,nfield,buried_x[k],buried_y[k],buried_z[k],cl,cl,cl*cl_factor,0)
            
        node_pos = np.append(node_pos,e_pt_idx) #add borehole nodes to electrode node positions 
        fh.write("Point{%s} In Surface{%i};//include node in meshing\n"%(str(e_pt_idx).strip('[').strip(']'), no_plane))
        fh.write("//end of buried electrodes.\n")
        elec_x_cache = np.append(elec_x_cache,buried_x)
        elec_z_cache = np.append(elec_z_cache,buried_z)
        # setFields(fh,nfield)
        
    if len(rem_idx)>0: #then we need to add remote electrodes to the mesh
        fh.write("\n//Remote electrodes \n")
        e_pt_idx = [0]*len(rem_idx)
        #put points in lower left hand corner of the outer region
        #if an electrode does have the remote flag, then chances are it has 
        #a bogus coordinate associated with it (ie. -99999)
        for k in range(len(rem_idx)):
            no_pts += 1
            remote_x = x_pts[0]-flank + (cl*cl_factor*(k+1))
            remote_z = b_max_depth  + (cl*cl_factor)
            fh.write("Point(%i) = {%.2f,%.2f,%.2f,cln};//remote electrode\n"%(no_pts,remote_x,0,remote_z))
            original_x[rem_idx[k]] = remote_x 
            original_z[rem_idx[k]] = remote_z 
            e_pt_idx[k] = no_pts
            elec_x_cache = np.append(elec_x_cache,original_x[rem_idx[k]])
            elec_z_cache = np.append(elec_z_cache,original_z[rem_idx[k]])
            
        node_pos = np.append(node_pos,e_pt_idx) #add remote electrode nodes to electrode node positions 
        fh.write("Point{%s} In Surface{1};\n"%(str(e_pt_idx).strip('[').strip(']')))
        fh.write('//End of remote electrodes.\n')
        
    if 'refine' in geom_input.keys():
        fh.write('\n//Refinement points \n')
        px = np.array(geom_input['refine'][0])
        pz = np.array(geom_input['refine'][1])
        if len(geom_input['refine']) == 3: 
            pl = np.array(geom_input['refine'][2])
        else:
            pl = np.ones(px.shape[0])*cl
        #work out which nodes are inside the fine mesh region 
        fmrx = np.append(x_pts,np.flipud(x_base))
        fmrz = np.append(z_pts,np.flipud(z_base))
        path = mpath.Path(np.c_[fmrx,fmrz])
        inside = path.contains_points(np.c_[px,pz])
        r_pt_idx = []
        for k in range(len(px)):
            if not inside[k]:
                continue 
            no_pts += 1 
            fh.write("Point(%i) = {%.2f,%.2f,%.2f,%.2f};//refinement point %i\n"%(no_pts, px[k], 0, pz[k], pl[k], k+1))
            r_pt_idx.append(no_pts) 
        fh.write("Point{%s} In Surface{%i};//include nodes in meshing\n"%(str(r_pt_idx).strip('[').strip(']'), no_plane))
        fh.write('//End refinement points\n')
                    
    fh.write("\n//End gmsh script\n")
    fh.close()
    dump("writing .geo to file completed, save location:\n%s\n"%os.getcwd())
    
    if len(node_pos) != len(elec_x_cache):
        warnings.warn("looks like something has gone wrong with node orderings, total x != total nodes.")
    
    #sort node ordering back into the original input ordering    
    #find the original indexes  
    original_idx = [0]*len(node_pos)
    for i in range(len(node_pos)):
        idx = (np.abs(elec_x_cache - original_x[i])<1e-15) & (np.abs(elec_z_cache - original_z[i])<1e-15)
        idx = idx.tolist()
        original_idx[i] = idx.index(True)

    ordered_node_pos = node_pos[original_idx].astype(int)
    
    return ordered_node_pos 
    
    
#%% 3D half space 
def halfspace3d(elec_x, elec_y, elec_z = None,
                fmd=-1, file_path='mesh3d.geo',
                cl=-1, cl_factor=-1, cln_factor=1e6, dp_len=-1, 
                mesh_refinement=None, use_fields=True, dump=None,
                flank_fac = 12):
    """
    Writes a gmsh .geo for a 3D half space with no topography for a hemispherical 
    mesh. Ignores the type of electrode. Z coordinates should be given as depth 
    below the surface! If Z != 0 then its assumed that the electrode is buried. 
    
    Parameters
    ----------
    elec_x: array like
        electrode x coordinates 
    elec_y: array like 
        electrode y coordinates 
    elec_z: array like 
        electrode z coordinates 
    fmd: float, optional 
        Depth of investigation of the survey. 
    file_path: string, optional 
        name of the generated gmsh file (can include file path also) (optional)
    cl: float, optional
        characteristic length (optional) of electrodes, essentially describes how big the nodes 
        assocaited elements will be on the electrodes. Usually no bigger than 5. If set as -1 (default)
        a characteristic length 1/4 the minimum electrode spacing is computed.
    cl_factor: float, optional 
        Characteristic length mulitplier for the sub-surface points on away from the 
        electrodes on the top of the fine mesh region.  This allows for tuning of 
        the incremental size increase with depth in the mesh, usually set to 2 such 
        that the elements at the DOI are twice as big as those
        at the surface. The reasoning for this is because the sensitivity of ERT drops
        off with depth. 
    cln_factor: float, optional
        Characteristic length mulitplier for the Neumann boundary points in the 
        coarse mesh region. 
    mesh_refinement: dict, pd.DataFrame, optional 
        Coordinates for discrete points in the mesh (advanced use cases). 
    dump : function, optional
        If None, output is printed using `print()`. Else the given function is passed.
    flank_factor: float, int, optional
        Defines how far away the outer radius of the of the mesh. The radius is 
        flank_factor * dp_len. 
    
    Returns
    ----------
    Node_pos: numpy array
        The indexes for the mesh nodes corresponding to the electrodes input, the ordering of the nodes
        should be the same as the input of the electrodes. 
    """
    
    if dump is None:
        def dump(x):
            print(x)
    if dp_len != -1 and dp_len<0:
        raise ValueError('Can not have a negative dipole length')
    if fmd != -1 and fmd<0:#then set to a default 
        raise ValueError('Can not have a negative depth of fine mesh')
    if cl != -1 and cl<0:
        raise ValueError('Can not have a negative mesh characteristic length')
        
    if file_path.find('.geo')==-1:
        file_path=file_path+'.geo'#add file extension if not specified already
    
    # check for z coordinates 
    if elec_z is None:
        elec_z = [0]*len(elec_x)

    # definition variable setup  
    nelec = len(elec_x)
    no_pts = 0
    nfield = 0 
    node_pos = np.zeros(nelec,dtype=int)
    dist = np.unique(find_dist(elec_x, elec_y, elec_z)) 
    
    dmin = dist[1]
    dmax = dist[-1]
    xmean = np.mean(elec_x)
    ymean = np.mean(elec_y)
    
    if cl==-1:
        cl = dmin/4 # characteristic length is 1/4 the minimum electrode distance

    if cl_factor == -1: 
        cl_factor = 5 
        
    if dp_len == -1: # compute largest possible dipole length 
        dp_len = dmax # maximum possible dipole length
    
    if fmd == -1: # compute depth of investigation
        fmd = dmax/3 # maximum possible dipole length / 3

    if fmd < abs(np.min(elec_z)):
        warnings.warn('depth of fine mesh is shallower than lowest electrode, adjusting...')
        fmd = abs(np.min(elec_z)) + (dp_len*2)

    outer_rad = dp_len*flank_fac # outer raduis of hemi-sphere 
    
    # create file 
    fh = open(file_path,'w')
    fh.write('//3D Half space (tetra) mesh for ResIPy\n')
    fh.write("Mesh.Binary = 0;//specify we want ASCII format\n")
    fh.write('SetFactory("OpenCASCADE");\n')
    fh.write('cl=%f;//electrode characteristic length\n'%float(cl))
    fh.write('cl_factor=%f;//Fine region characteristic factor\n'%float(cl_factor))
    fh.write('cln_factor=%f;//Neumann characteristic factor\n'%float(cln_factor))
    fh.write('//Create a hemisphere which approximates an infinite halfspace\n')
    
    template = 'Sphere(1) = {%f, %f, %f, %f, -Pi/2, 0, Pi*2};//create a sphere\n' 
    no_pts+=2 # seems adding a sphere adds in 2 points 
    fh.write(template%(xmean,ymean,0,outer_rad))#
    fh.write("Mesh.CharacteristicLengthFromPoints=1;\n")
    # fh.write("MeshSize{1,2} = cl*cln_factor;")
    
    # create electrodes
    fh.write('//Start electrodes, include refinement sphere around the electrode\n')
    template = "Point(%i) = {%f, %f, %f, cl};//electrode coordinate\n"
    for i in range(nelec):
        no_pts += 1 
        fh.write(template%(no_pts,elec_x[i],elec_y[i],elec_z[i]))
        node_pos[i] = no_pts
        if elec_z[i] == 0:
            fh.write("Point{%i} In Surface{2};//specify electrode surface\n"%no_pts)
        else:
            fh.write("Point{%i} In Volume{1};//specify electrode volume\n"%no_pts)
        if use_fields: 
            nfield += 1 
            addBallField(fh,nfield,elec_x[i],elec_y[i],elec_z[i],
                         dmin/4,0)
    if use_fields:
        setFields(fh, nfield)
    fh.write("//End of electrodes\n")
    
    

    if mesh_refinement is None and min(elec_z) == 0: 
        
        # add some safety net and add some refinement at depth for the user 
        # if no mesh refinement provided and only surface electrodes are present 
        
        ref_x = [min(elec_x), min(elec_x), max(elec_x), max(elec_x), xmean]
        ref_y = [min(elec_y), max(elec_y), min(elec_y), max(elec_y), ymean]
        tree = cKDTree(np.c_[ref_x,ref_y])
        idist, _ = tree.query(np.c_[ref_x,ref_y],2)
        tokeep = [True]*5 
        for i in range(5): 
            if min(idist[i,1]) < 1e-16: 
                tokeep[i] = False 
        ref_x = [ref_x[b] for b in tokeep]
        ref_y = [ref_y[b] for b in tokeep]

        fh.write("//Subsurface refinement fields\n")
        template = "Point(%i) = {%f, %f, %f, cl*cl_factor};//refinement coordinate\n"
        for i in range(len(ref_x)):
            no_pts += 1 
            fh.write(template%(no_pts,ref_x[i],ref_y[i],-fmd))
            fh.write("Point{%i} In Volume{1};//specify refinement in volume\n"%no_pts)
        
        
    # check if any mesh refinement is requested 
    elif mesh_refinement is not None:
        fh.write('//Custom refinement points\n')
        # find surface points 
        rx = mesh_refinement['x']
        ry = mesh_refinement['y']
        rz = mesh_refinement['z']
        npoints = len(rx) # number of refinement points 
        if 'cl' in mesh_refinement.keys():
            rcl = mesh_refinement['cl'] # use cl specified for each point 
            fh.write('//Refinement point characteristic length specified for each point individually\n')
        else:
            rcl = [cl]*npoints # use same cl as for electrodes 
            fh.write('//Refinement point characteristic length the same as for electrodes\n')
        for i in range(npoints):
            no_pts += 1
                
            if rz[i] == 0: # it is on the surface! 
                fh.write("Point (%i) = {%.16f, %.16f, %.16f, %.16f};\n"%(no_pts, rx[i], ry[i], 0, rcl[i]))
                fh.write("Point{%i} In Surface{%i};//surface refinement point\n"%(no_pts,2))# put the point on surface
            else:
                if rz[i]>0:
                    fh.close()
                    raise ValueError("point z coordinate is greater than 0 in gmshWrap.py and you can't have refinement points above the surface!")
                fh.write("Point (%i) = {%.16f, %.16f, %.16f, %.16f};\n"%(no_pts, rx[i], ry[i], rz[i], rcl[i]))
                fh.write("Point{%i} In Volume{%i};//buried refinement point\n"%(no_pts,1))# put the point in volume 
        fh.write('//End mesh refinement points\n')
        
    fh.write("//End of gmsh script\n")
        
    fh.close()
    
    return node_pos

def __cylinder25d(elec_x, elec_y, elec_z = None,
                  fmd=-1, file_path='mesh3d.geo',
                  cl=-1, cl_factor=-1, cln_factor=1000, dp_len=-1, 
                  mesh_refinement=None, use_fields=False, dump=None,
                  flank_fac = 10):
    """
    Create a cylindrical half space in the case of long 2D lines. 

    """

    if file_path.find('.geo')==-1:
        file_path=file_path+'.geo'#add file extension if not specified already
    
    nelec = len(elec_x)
        
    no_pts = 0
    nfield = 0 
    no_lines = 0 
    no_surf = 0 
    node_pos = np.zeros(nelec,dtype=int)
    dist = np.unique(find_dist(elec_x, np.zeros(nelec), elec_z)) 
    
    dmin = dist[1]
    dmax = dist[-1]

    outer_rad = dmax*flank_fac
    
    # create file 
    fh = open(file_path,'w')
    fh.write('//3D Half space (tetra) mesh for ResIPy which is a 2D line\n')
    fh.write("Mesh.Binary = 0;//specify we want ASCII format\n")
    fh.write('SetFactory("OpenCASCADE");\n')
    fh.write('cl=%f;//electrode characteristic length\n'%float(cl))
    fh.write('cl_factor=%f;//Fine region characteristic factor\n'%float(cl_factor))
    fh.write('cln_factor=%f;//Neumann characteristic factor\n'%float(cln_factor))
    fh.write('r=%f;\n'%float(outer_rad))
    
    # create electrodes
    fh.write('//Start electrodes, include refinement sphere around the electrode\n')
    template = "Point(%i) = {%f, %f, %f, cl};//electrode\n"
    pts_cache = []
    lns_cache = []
    bur_cache = [] # buried electrodes 
    sur_cache = [] # surface cache 
    for i in range(len(elec_x)):
        no_pts += 1 
        fh.write(template%(no_pts,elec_x[i],elec_y[0],elec_z[i]))
        node_pos[i] = no_pts
        if elec_z[i] == 0:
            pts_cache.append(no_pts)
        else:
            bur_cache.append(no_pts)
        nfield += 1 
        addBallField(fh,nfield,elec_x[i],elec_y[0],elec_z[i],
                     dmin/4)
    setFields(fh, nfield)
    
    for i in range(1,len(pts_cache)):
        no_lines += 1 
        fh.write("Line(%i) = {%i, %i};\n"%(no_lines,pts_cache[i-1],pts_cache[i]))
        lns_cache.append(no_lines)
        
    template = "Point(%i) = {%f, %f, %f, cl*cl_factor};//depth refinement point\n"
    no_pts += 1 
    fh.write(template%(no_pts,elec_x[0],elec_y[0],-fmd))
    no_pts += 1 
    fh.write(template%(no_pts,elec_x[-1],elec_y[0],-fmd))
    
    no_lines += 1 
    fh.write("Line(%i) = {%i, %i};\n"%(no_lines,pts_cache[-1],no_pts))
    lns_cache.append(no_lines)
    
    no_lines += 1 
    fh.write("Line(%i) = {%i, %i};\n"%(no_lines,no_pts,no_pts-1))
    lns_cache.append(no_lines)
    
    no_lines += 1 
    fh.write("Line(%i) = {%i, %i};\n"%(no_lines,pts_cache[0],no_pts-1))
    lns_cache.append(no_lines)
    
    no_surf +=1 
    fh.write('Curve Loop(%i) = {'%no_surf)
    for i in range(len(lns_cache)-1):
        fh.write('%i, '%lns_cache[i])
    fh.write('-%i};\n'%lns_cache[-1])
    
    fh.write('Plane Surface(%i) = {%i};\n'%(no_surf,no_surf))

    # specify points below surface are in 2D mesh region         
    for pt in bur_cache: 
        fh.write('Point{%i} In Surface{1};\n'%pt)
        
    # add lines towards the right flanks of the survey 
    right_pts = [0]*4 
    template = "Point(%i) = {%f+r, %f, %f, cl*cln_factor};//outer boundary point\n"
    no_pts += 1 
    fh.write(template%(no_pts,elec_x[-1],elec_y[0],0))
    right_pts[0] = no_pts 
    
    template = "Point(%i) = {%f+r, %f-r, %f, cl*cln_factor};//outer boundary point\n"
    no_pts += 1 
    fh.write(template%(no_pts,elec_x[-1],elec_y[0],0))
    right_pts[1] = no_pts 
    
    template = "Point(%i) = {%f+r, %f+r, %f, cl*cln_factor};//outer boundary point\n"
    no_pts += 1 
    fh.write(template%(no_pts,elec_x[-1],elec_y[0],0))
    right_pts[2] = no_pts 
    
    template = "Point(%i) = {%f+r, %f, -r, cl*cln_factor};//outer boundary point\n"
    no_pts += 1 
    fh.write(template%(no_pts,elec_x[-1],elec_y[0]))
    right_pts[3] = no_pts 
    
    # add lines towards the left flanks of the survey 
    left_pts = [0]*4 
    template = "Point(%i) = {%f-r, %f, %f, cl*cln_factor};//outer boundary point\n"
    no_pts += 1 
    fh.write(template%(no_pts,elec_x[0],elec_y[0],0))
    left_pts[0] = no_pts 
    
    template = "Point(%i) = {%f-r, %f-r, %f, cl*cln_factor};//outer boundary point\n"
    no_pts += 1 
    fh.write(template%(no_pts,elec_x[0],elec_y[0],0))
    left_pts[1] = no_pts 
    
    template = "Point(%i) = {%f-r, %f+r, %f, cl*cln_factor};//outer boundary point\n"
    no_pts += 1 
    fh.write(template%(no_pts,elec_x[0],elec_y[0],0))
    left_pts[2] = no_pts 

    template = "Point(%i) = {%f-r, %f, -r, cl*cln_factor};//outer boundary point\n"
    no_pts += 1 
    fh.write(template%(no_pts,elec_x[0],elec_y[0]))
    left_pts[3] = no_pts 
    
    
    # start to construct the surfaces of the cylinder 
    fh.write("//Following lines construct the shell of the halfspace cylinder\n")
    curve_lines_0 = [0]*4 
    curve_lines_1 = [0]*4 
    top_surf_line = [0]*6 
    right_surf_ln = [0]*4 
    left_surf_lin = [0]*4 
    
    # add circles 
    template = "Circle(%i) = {%i, %i, %i};\n"
    no_lines +=1 
    fh.write(template%(no_lines,right_pts[1],right_pts[0],right_pts[3]))
    curve_lines_0[0] = no_lines 
    right_surf_ln[2] = no_lines 
    
    no_lines +=1  
    fh.write(template%(no_lines,right_pts[2],right_pts[0],right_pts[3]))
    curve_lines_1[0] = no_lines 
    right_surf_ln[3] = no_lines 
    
    no_lines +=1 
    fh.write(template%(no_lines,left_pts[1],left_pts[0],left_pts[3]))
    curve_lines_0[2] = no_lines 
    left_surf_lin[2] = no_lines 

    no_lines +=1  
    fh.write(template%(no_lines,left_pts[2],left_pts[0],left_pts[3]))
    curve_lines_1[2] = no_lines 
    left_surf_lin[3] = no_lines 
    
    
    # add more lines to form top surface 
    template = "Line(%i) = {%i, %i};\n"
    no_lines +=1  
    fh.write(template%(no_lines,right_pts[2],left_pts[2]))
    curve_lines_1[1] = no_lines
    top_surf_line[0] = no_lines 

    no_lines +=1  
    fh.write(template%(no_lines,right_pts[1],left_pts[1]))
    curve_lines_0[1] = no_lines
    top_surf_line[3] = no_lines 

    no_lines +=1  
    fh.write(template%(no_lines,right_pts[0],right_pts[1]))
    right_surf_ln[1] = no_lines 
    top_surf_line[2] = no_lines 

    no_lines +=1  
    fh.write(template%(no_lines,right_pts[0],right_pts[2]))
    right_surf_ln[0] = no_lines 
    top_surf_line[1] = no_lines 
    
    no_lines +=1  
    fh.write(template%(no_lines,left_pts[0],left_pts[1]))
    left_surf_lin[1] = no_lines 
    top_surf_line[4] = no_lines 

    no_lines +=1  
    fh.write(template%(no_lines,left_pts[0],left_pts[2]))
    left_surf_lin[0] = no_lines 
    top_surf_line[5] = no_lines 
    
    # add line at base of cylinder 
    template = "Line(%i) = {%i, %i};\n"
    no_lines +=1  
    fh.write(template%(no_lines,right_pts[3],left_pts[3]))
    curve_lines_0[3] = no_lines 
    curve_lines_1[3] = no_lines 
    
    # add the curved surfaces 
    template = "Curve Loop(%i) = {%i, -%i, -%i, %i};//cylinder half space surface\n"
    no_surf += 1 
    fh.write(template%(no_surf,curve_lines_0[0],curve_lines_0[1],curve_lines_0[2],curve_lines_0[3]))
    fh.write("Surface(%i) = {%i};\n"%(no_surf,no_surf))
    sur_cache.append(no_surf)
    no_surf += 2 # i dont know why but open cascade seems to add 2 surfaces in the background per curve loop 
    fh.write(template%(no_surf,curve_lines_1[0],curve_lines_1[1],curve_lines_1[2],curve_lines_1[3]))
    fh.write("Surface(%i) = {%i};\n"%(no_surf,no_surf))
    sur_cache.append(no_surf)
    
    # add end on surfaces 
    template = "Curve Loop(%i) = {%i, %i, %i, %i};//cylinder end on surface\n"
    no_surf += 2 
    fh.write(template%(no_surf,right_surf_ln[0],right_surf_ln[1],right_surf_ln[2],right_surf_ln[3]))
    fh.write("Plane Surface(%i) = {%i};\n"%(no_surf,no_surf))
    sur_cache.append(no_surf)
    no_surf += 2 
    fh.write(template%(no_surf,left_surf_lin[0],left_surf_lin[1],left_surf_lin[2],left_surf_lin[3]))
    fh.write("Plane Surface(%i) = {%i};\n"%(no_surf,no_surf))
    sur_cache.append(no_surf)
    
    #add top surface 
    no_surf += 2 
    fh.write("Curve Loop(%i) = {"%no_surf)
    for i in range(5):
        fh.write("%i ,"%top_surf_line[i])
    fh.write("%i};\n"%top_surf_line[-1])
    fh.write("Plane Surface(%i) = {%i};\n"%(no_surf,no_surf))
    sur_cache.append(no_surf)
    
    # put the electrodes into the top surface 
    fh.write("//Add the lines of the 2D survey to the top surface of the mesh\n")
    for i in range(len(lns_cache)-3):
        line = lns_cache[i]
        fh.write("Line{%i} In Surface {%i};\n"%(line,no_surf))
        
    # make volume 
    fh.write("//Combine surfaces together to make the meshing volume\n")
    fh.write("Surface Loop(1) = {")
    for i in range(len(sur_cache)-1):
        fh.write("%i ,"%sur_cache[i])
    fh.write("%i};\n"%sur_cache[-1])
    fh.write("Volume(1) = {1};\n")
    
    # finally add the 2d surface into the volume 
    fh.write("Surface{1} In Volume{1};\n")
    
    fh.close()
    
    return node_pos 
    
def halfspace3dline2d(elec_x, elec_y, elec_z = None,
                      fmd=-1, file_path='mesh3d.geo',
                      cl=-1, cl_factor=-1, cln_factor=1000, dp_len=-1, 
                      mesh_refinement=None, use_fields=False, dump=None,
                      flank_fac = 12, geom='hemisphere'):
    """
    Writes a gmsh .geo for a 3D half space with no topography for a 2D line in a 
    3D half space problem. Either a hemispherical or cylindrical mesh can be used. 
    mesh for a survey setup that is a single 2D line. Ignores the type of electrode. 
    Z coordinates should be given as depth below the surface! If Z != 0 then 
    its assumed that the electrode is buried. 
    
    Parameters
    ----------
    elec_x: array like
        electrode x coordinates 
    elec_z: array like 
        electrode z coordinates 
    fmd: float, optional 
        Depth of investigation of the survey. 
    file_path: string, optional 
        name of the generated gmsh file (can include file path also) (optional)
    cl: float, optional
        characteristic length (optional) of electrodes, essentially describes how big the nodes 
        assocaited elements will be on the electrodes. Usually no bigger than 5. If set as -1 (default)
        a characteristic length 1/4 the minimum electrode spacing is computed.
    cl_factor: float, optional 
        Characteristic length mulitplier for the sub-surface points on away from the 
        electrodes on the top of the fine mesh region.  This allows for tuning of 
        the incremental size increase with depth in the mesh, usually set to 2 such 
        that the elements at the DOI are twice as big as those
        at the surface. The reasoning for this is because the sensitivity of ERT drops
        off with depth. 
    cln_factor: float, optional
        Characteristic length mulitplier for the Neumann boundary points in the 
        coarse mesh region. 
    mesh_refinement: dict, pd.DataFrame, optional 
        Coordinates for discrete points in the mesh (advanced use cases). 
    dump : function, optional
        If None, output is printed using `print()`. Else the given function is passed.
    flank_fac: float, int, optional
        Defines how far away the outer radius of the of the mesh. The radius is 
        flank_factor * dp_len. 
    geom: str, optional 
        Change from hemisphere to make a cylindrical mesh!
    
    Returns
    ----------
    Node_pos: numpy array
        The indexes for the mesh nodes corresponding to the electrodes input, the ordering of the nodes
        should be the same as the input of the electrodes. 
    """
    if dump is None:
        def dump(x):
            print(x)
    if dp_len != -1 and dp_len<0:
        raise ValueError('Can not have a negative dipole length')
    if fmd != -1 and fmd<0:#then set to a default 
        raise ValueError('Can not have a negative depth of fine mesh')
    if cl != -1 and cl<0:
        raise ValueError('Can not have a negative mesh characteristic length')
    
    if file_path.find('.geo')==-1:
        file_path=file_path+'.geo'#add file extension if not specified already
    
    nelec = len(elec_x)
    if elec_z is None:
        elec_z = [0]*len(elec_x)
        
    no_pts = 0
    nfield = 0 
    no_lines = 0 
    node_pos = np.zeros(nelec,dtype=int)
    dist = np.unique(find_dist(elec_x, np.zeros(nelec), elec_z)) 
    
    dmin = dist[1]
    dmax = dist[-1]
    xmean = np.mean(elec_x)
    ymean = elec_y[0]

    if cl==-1:
        cl = dmin/4 # characteristic length is 1/4 the minimum electrode distance
    
    if cl_factor == -1: 
        cl_factor = 5 
        
    if dp_len == -1: # compute largest possible dipole length 
        dp_len = dmax # maximum possible dipole length
    
    if fmd == -1: # compute depth of investigation
        fmd = dmax/3 # maximum possible dipole length / 3
    
    if fmd < abs(np.min(elec_z)):
        warnings.warn('depth of fine mesh is shallower than lowest electrode, adjusting...')
        fmd = abs(np.min(elec_z)) + (dp_len*2)
        
    if geom != 'hemisphere':
        node_pos = __cylinder25d(elec_x, elec_y, elec_z, fmd, file_path,cl,
                                 cl_factor,cln_factor,dp_len,mesh_refinement,
                                 use_fields,dump, flank_fac)
        return node_pos 

    outer_rad = dmax*flank_fac
    
    # create file 
    fh = open(file_path,'w')
    fh.write('//3D Half space (tetra) mesh for ResIPy which is a 2D line\n')
    fh.write("Mesh.Binary = 0;//specify we want ASCII format\n")
    fh.write('SetFactory("OpenCASCADE");\n')
    fh.write('cl=%f;//electrode characteristic length\n'%float(cl))
    fh.write('cl_factor=%f;//Fine region characteristic factor\n'%float(cl_factor))
    fh.write('cln_factor=%f;//Neumann characteristic factor\n'%float(cln_factor))
    
    template = 'Sphere(1) = {%f, %f, %f, %f, -Pi/2, 0, Pi*2};//create a sphere\n' 
    no_pts+=2 # seems adding a sphere adds in 2 points 
    no_lines+=3 # and three lines 
    fh.write(template%(xmean,ymean,0.0,outer_rad))
    
    # create electrodes
    fh.write('//Start electrodes, include refinement sphere around the electrode\n')
    template = "Point(%i) = {%f, %f, %f, cl};//electrode\n"
    pts_cache = []
    lns_cache = []
    bur_cache = [] # buried electrodes 
    for i in range(len(elec_x)):
        no_pts += 1 
        fh.write(template%(no_pts,elec_x[i],elec_y[0],elec_z[i]))
        node_pos[i] = no_pts
        if elec_z[i] == 0:
            fh.write("Point{%i} In Surface{2};\n"%no_pts)
            pts_cache.append(no_pts)
        else:
            # fh.write("Point{%i} In Volume{1};\n"%no_pts)
            bur_cache.append(no_pts)
        nfield += 1 
        addBallField(fh,nfield,elec_x[i],elec_y[0],elec_z[i],
                     dmin/4)
    setFields(fh, nfield)
    
    for i in range(1,len(pts_cache)):
        no_lines += 1 
        fh.write("Line(%i) = {%i, %i};\n"%(no_lines,pts_cache[i-1],pts_cache[i]))
        lns_cache.append(no_lines)
        
    template = "Point(%i) = {%f, %f, %f, cl*cl_factor};//depth refinement point\n"
    no_pts += 1 
    fh.write(template%(no_pts,elec_x[0],elec_y[0],-fmd))
    no_pts += 1 
    fh.write(template%(no_pts,elec_x[-1],elec_y[0],-fmd))
    
    no_lines += 1 
    fh.write("Line(%i) = {%i, %i};\n"%(no_lines,pts_cache[-1],no_pts))
    lns_cache.append(no_lines)
    
    no_lines += 1 
    fh.write("Line(%i) = {%i, %i};\n"%(no_lines,no_pts,no_pts-1))
    lns_cache.append(no_lines)
    
    no_lines += 1 
    fh.write("Line(%i) = {%i, %i};\n"%(no_lines,pts_cache[0],no_pts-1))
    lns_cache.append(no_lines)
    
    fh.write('Curve Loop(3) = {')
    for i in range(len(lns_cache)-1):
        fh.write('%i, '%lns_cache[i])
    fh.write('-%i};\n'%lns_cache[-1])
    
    fh.write('Plane Surface(3) = {3};\n')
    fh.write('Surface{3} In Volume{1};\n')
    
    for i in range(len(lns_cache)-3):
        fh.write('Line{%i} In Surface{2};\n'%lns_cache[i])

    # specify points below surface are in 2D mesh region         
    for pt in bur_cache: 
        fh.write('Point{%i} In Surface{3};\n'%pt)
        
    fh.close()
    
    return node_pos 

def box_3d(electrodes, padding=20, fmd=-1, file_path='mesh3d.geo',
           cl=-1, cl_corner=-1, cl_factor=8, cln_factor=100, dp_len=-1, 
           mesh_refinement=None, use_fields=False, dump=None, 
           ball_refinement=None, add_veroni_refinement=None):
    """
    writes a gmsh .geo for a 3D half space with no topography. Ignores the type of electrode. 
    Z coordinates should be given as depth below the surface! If Z != 0 then its assumed that the
    electrode is buried. 
    
    !! Depreciated code as of Dec 2023 !! 
    
    Parameters
    ----------
    electrodes: list of array likes
        first column/list is the x coordinates of electrode positions, second column
        is the elevation. Z coordinates must normalised to the flat surface if given.ie. 
        z is the depth below the surface. 
    padding: float, optional
        Padding in percent on the size the fine mesh region extent. Must be bigger than 0.
    fmd: float, optional 
        Depth of investigation of the survey. 
    file_path: string, optional 
        name of the generated gmsh file (can include file path also) (optional)
    cl: float, optional
        characteristic length (optional) of electrodes, essentially describes how big the nodes 
        assocaited elements will be on the electrodes. Usually no bigger than 5. If set as -1 (default)
        a characteristic length 1/4 the minimum electrode spacing is computed.
    cl_corner: float, optional
        Characteristic length mulitplier for the surface points on away from the 
        electrodes on the top of the fine mesh region. The reasoning for this 
        is because the sensitivity of ERT drops off with distance from electrodes. 
    cl_factor: float, optional 
        Characteristic length mulitplier for the sub-surface points on away from the 
        electrodes on the top of the fine mesh region.  This allows for tuning of 
        the incremental size increase with depth in the mesh, usually set to 2 such 
        that the elements at the DOI are twice as big as those
        at the surface. The reasoning for this is because the sensitivity of ERT drops
        off with depth. 
    cln_factor: float, optional
        Characteristic length mulitplier for the nuemmon boundary points in the 
        coarse mesh region. 
    mesh_refinement: dict, pd.DataFrame, optional 
        Coordinates for discrete points in the mesh (advanced use cases). 
    dump : function, optional
        If None, output is printed using `print()`. Else the given function is passed.
    
    Returns
    ----------
    Node_pos: numpy array
        The indexes for the mesh nodes corresponding to the electrodes input, the ordering of the nodes
        should be the same as the input of 'electrodes'
    .geo: file
        Can be run in gmsh

    NOTES
    ----------            
    electrodes format: 
        
            electrodes = [[x1,x2,x3,...],[y1,y2,y3,...],[z1,z2,z3,...]]
                
    #### TODO: search through each set of points and check for repeats ?
    """
    if dump is None:
        def dump(x):
            print(x)
    if dp_len != -1 and dp_len<0:
        raise ValueError('Can not have a negative dipole length')
    if fmd != -1 and fmd<0:#then set to a default 
        raise ValueError('Can not have a negative depth of fine mesh')
    if cl != -1 and cl<0:
        raise ValueError('Can not have a negative mesh characteristic length')
    
    elec_x = electrodes[0]
    elec_y = electrodes[1]
    try: 
        elec_z = electrodes[2]
    except IndexError:
        elec_z = [0]*len(elec_x)
    
    if len(elec_x) != len(elec_y):
        raise ValueError("The length of the x coordinate array does not match of the Y coordinate")
        
    if file_path.find('.geo')==-1:
        file_path=file_path+'.geo'#add file extension if not specified already
    
    dist_sort = np.unique(find_dist(elec_x,elec_y,elec_z))
    if cl==-1:
        cl = dist_sort[1]/2 # characteristic length is 1/2 the minimum electrode distance

    if cl_corner == -1: 
        cl_corner = 3*cl 
        
    if dp_len == -1: # compute largest possible dipole length 
        dp_len = dist_sort[-1] # maximum possible dipole length
    
    if fmd == -1: # compute depth of investigation
        fmd = dist_sort[-1]/3 # maximum possible dipole length / 3

    if fmd < abs(np.min(elec_z)):
        warnings.warn('depth of fine mesh is shallower than lowest electrode, adjusting...')
        fmd = abs(np.min(elec_z)) + (dp_len*2)
        
    # print('fmd in gmshWrap.py: %f'%fmd)
    ### start to write to file ### 
    fh = open(file_path,'w') #file handle
    
    fh.write("//3D half space problem mesh for ResIPy - no topography\n")
    fh.write("Mesh.Binary = 0;//specify we want ASCII format\n")
    fh.write("cl=%.2f;//define characteristic length for fine mesh region\n" %cl)
    
    #create square around all of the electrodes
    x_dist = abs(np.max(elec_x) - np.min(elec_x))
    y_dist = abs(np.max(elec_y) - np.min(elec_y))
    if x_dist<0.2:x_dist=5 # protection against small (or zero) padding 
    if y_dist<0.2:y_dist=5
        
    max_x = np.max(elec_x) + (padding/100)*x_dist
    min_x = np.min(elec_x) - (padding/100)*x_dist
    max_y = np.max(elec_y) + (padding/100)*y_dist
    min_y = np.min(elec_y) - (padding/100)*y_dist
    
    fh.write("//Fine mesh region.\n")
    #add points to file at z = zero 
    no_pts = 1
    loop_pt_idx=[no_pts]
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, %.2f};\n"%(no_pts, max_x, max_y, 0, cl*cl_corner))#multiply by 2 grow the elements away from electrodes 
    no_pts += 1
    loop_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, %.2f};\n"%(no_pts, max_x, min_y, 0, cl*cl_corner))
    no_pts += 1
    loop_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, %.2f};\n"%(no_pts, min_x, min_y, 0, cl*cl_corner))
    no_pts += 1
    loop_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, %.2f};\n"%(no_pts, min_x, max_y, 0, cl*cl_corner))
    
    #add line loop
    no_lns = 0 
    for i in range(4):
        no_lns += 1 
        if i == 3:
            fh.write("Line(%i) = {%i,%i};\n"%(no_lns,loop_pt_idx[i],loop_pt_idx[0]))
        else:
            fh.write("Line(%i) = {%i,%i};\n"%(no_lns,loop_pt_idx[i],loop_pt_idx[i+1]))
    
    fmd = abs(fmd)
    #add points below surface to make a rectangular volume  
    fh.write("cl2=cl*%.2f;//define characteristic length for base of fine mesh region\n" %cl_factor)       
    no_pts += 1
    loop2_pt_idx=[no_pts]
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl2};\n"%(no_pts, max_x, max_y, -fmd))
    no_pts += 1
    loop2_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl2};\n"%(no_pts, max_x, min_y, -fmd))
    no_pts += 1
    loop2_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl2};\n"%(no_pts, min_x, min_y, -fmd))
    no_pts += 1
    loop2_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl2};\n"%(no_pts, min_x, max_y, -fmd))
    
    #add line loops to connect the points just made 
    for i in range(4):
        no_lns += 1 
        if i == 3:
            fh.write("Line(%i) = {%i,%i};\n"%(no_lns,loop2_pt_idx[i],loop2_pt_idx[0]))
        else:
            fh.write("Line(%i) = {%i,%i};\n"%(no_lns,loop2_pt_idx[i],loop2_pt_idx[i+1]))

    #connect the top and bottom of the mesh 
    for i in range(4):
        no_lns += 1 
        fh.write("Line(%i) = {%i,%i};\n"%(no_lns,loop_pt_idx[i],loop2_pt_idx[i]))   
        
    fh.write("//End fine mesh region points\n" )
        
    #Neumann boundary 
    flank_x = 5 * dp_len
    flank_y = 5 * dp_len
    if any(elec_z < 0):       
        flank_z = 5 * dp_len
        fh.write('//Maximum mesh depth calculated as 5*max dipole length (%f)\n'%dp_len)
    else:
        flank_z = fmd + (3*fmd)
        fh.write('//Maximum mesh depth calculated as fmd  + (3*fmd) (fmd=%f)\n'%fmd)
        
    fh.write("//Neumannn boundary points\n")
    cln = cl*cln_factor # nuemom boundary characteristic length 
    fh.write("cln = %.2f;//characteristic length for background region\n"%cln)
    no_pts += 1
    nmn_pt_idx=[no_pts]
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cln};\n"%(no_pts, max_x+flank_x, max_y+flank_y, 0))
    no_pts += 1
    nmn_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cln};\n"%(no_pts, max_x+flank_x, min_y-flank_y, 0))
    no_pts += 1
    nmn_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cln};\n"%(no_pts, min_x-flank_x, min_y-flank_y, 0))
    no_pts += 1
    nmn_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cln};\n"%(no_pts, min_x-flank_x, max_y+flank_y, 0))
    
    for i in range(4):
        no_lns += 1 
        if i == 3:
            fh.write("Line(%i) = {%i,%i};\n"%(no_lns,nmn_pt_idx[i],nmn_pt_idx[0]))
        else:
            fh.write("Line(%i) = {%i,%i};\n"%(no_lns,nmn_pt_idx[i],nmn_pt_idx[i+1]))
    
    #base of background region      
    no_pts += 1
    nmn2_pt_idx=[no_pts]
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cln};\n"%(no_pts, max_x+flank_x, max_y+flank_y, -fmd - flank_z))
    no_pts += 1
    nmn2_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cln};\n"%(no_pts, max_x+flank_x, min_y-flank_y, -fmd - flank_z))
    no_pts += 1
    nmn2_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cln};\n"%(no_pts, min_x-flank_x, min_y-flank_y, -fmd - flank_z))
    no_pts += 1
    nmn2_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cln};\n"%(no_pts, min_x-flank_x, max_y+flank_y, -fmd - flank_z))
    
    for i in range(4):
        no_lns += 1 
        if i == 3:
            fh.write("Line(%i) = {%i,%i};\n"%(no_lns,nmn2_pt_idx[i],nmn2_pt_idx[0]))
        else:
            fh.write("Line(%i) = {%i,%i};\n"%(no_lns,nmn2_pt_idx[i],nmn2_pt_idx[i+1]))
            
    for i in range(4):
        no_lns += 1 
        fh.write("Line(%i) = {%i,%i};\n"%(no_lns,nmn_pt_idx[i],nmn2_pt_idx[i]))  
    
    fh.write("//End of nuemmon boundary points\n")
    #add relevant line loops for top surface      
    fh.write("Line Loop(1) = {1,2,3,4};\n") # fine mesh region top line loop         
    fh.write("Line Loop(2) = {13,14,15,16};\n")# top nuemmon boundary loop 
    fh.write("Plane Surface(1) = {1};\n")
    fh.write("Plane Surface(2) = {2,1};\n") # top mesh with a hole in it for fine mesh region 
    
    #add below ground fine mesh surfaces 
    fh.write("//Below ground fine mesh surfaces\n")
    fh.write("Line Loop(3) = {5,6,7,8};\n") # basal surface
    fh.write("Line Loop(4) = {5,-10,-1,9};\n") 
    fh.write("Line Loop(5) = {11,-6,-10,2};\n") 
    fh.write("Line Loop(6) = {11,7,-12,-3};\n") 
    fh.write("Line Loop(7) = {8,-9,-4,12};\n") 
    for i in range(3,8):
        fh.write("Plane Surface(%i) = {%i};\n"%(i,i))
    
    fh.write("Surface Loop (1) = {1,3,4,5,6,7};\n") # fine mesh region volume 
    fh.write("Volume(1) = {1};//End Fine mesh region surfaces.\n")
    
    fh.write("//Below ground background surfaces\n")
    fh.write("Line Loop(8) = {20, 17, 18, 19};\n")
    fh.write("Line Loop(9) = {24, -19, -23, 15};\n")
    fh.write("Line Loop(10) = {16, 21, -20, -24};\n")
    fh.write("Line Loop(11) = {17, -22, -13, 21};\n")
    fh.write("Line Loop(12) = {18, -23, -14, 22};\n")
    for i in range(8,13):
        fh.write("Plane Surface(%i) = {%i};\n"%(i,i))
        
    fh.write("Surface Loop (2) = {2,3,4,5,6,7,8,9,10,11,12};\n") # background mesh region volume 
    fh.write("Volume(2) = {2};//End background mesh surfaces\n")   
    
    # add electrodes to mesh 
    fh.write("//Start electrode positions.\n")
    node_pos = [0]*len(elec_x)
    for i in range(len(elec_x)):
        no_pts += 1
        node_pos[i] = no_pts
        if elec_z[i] == 0:
            fh.write("Point (%i) = {%.16f,%.16f,%.16f, cl};\n"%(no_pts, elec_x[i], elec_y[i], 0))
            fh.write("Point{%i} In Surface{1};//surface electrode\n"%(no_pts))# put the point on surface
        else:
            if elec_z[i]>0:#
                fh.close() # close file 
                raise ValueError("electrode z coordinate is greater than 0 in gmshWrap.py and you can't have buried electrodes above the surface!")
            fh.write("Point (%i) = {%.16f,%.16f,%.16f, cl};\n"%(no_pts, elec_x[i], elec_y[i], elec_z[i]))
            fh.write("Point{%i} In Volume{1};//buried electrode\n"%(no_pts))# put the point in volume 
    
    nfield=0 
    if use_fields: 
        fh.write('//Electrode refinement fields\n')
        for i in range(len(elec_x)):
            nfield += 1 
            add_ball_field(fh, nfield, elec_x[i], elec_y[i], elec_z[i],
                           dist_sort[1]/2, cl, cl*cl_corner,0)
    fh.write("//End electrodes\n")
    
    # check if any mesh refinement is requested 
    if mesh_refinement is not None:
        fh.write('//Start mesh refinement points\n')
        # find surface points 
        rx = mesh_refinement['x']
        ry = mesh_refinement['y']
        rz = mesh_refinement['z']
        npoints = len(rx) # number of refinement points 
        if 'cl' in mesh_refinement.keys():
            rcl = mesh_refinement['cl'] # use cl specified for each point 
            fh.write('//Refinement point characteristic length specified for each point individually\n')
        else:
            rcl = [cl]*npoints # use same cl as for electrodes 
            fh.write('//Refinement point characteristic length the same as for electrodes\n')
        for i in range(npoints):
            no_pts += 1
            # check if point in refined zone or not 
            if rx[i] > min_x and rx[i] < max_x and ry[i] > min_y and ry[i] < max_y and rz[i] > -fmd:
                sur_idx = 1
                vol_idx = 1 
            else: # otherwise point placed in mesh outside of fine mesh region (zone 1)
                continue 
                sur_idx = 2 
                vol_idx = 2 
            
            if rx[i] == min_x or rx[i] == max_x or ry[i] == min_y or ry[i] == max_y or rz[i] == -fmd:
                continue 
                
            if rz[i] == 0: # it is on the surface! 
                fh.write("Point (%i) = {%.16f, %.16f, %.16f, %.16f};\n"%(no_pts, rx[i], ry[i], 0, rcl[i]))
                fh.write("Point{%i} In Surface{%i};//surface refinement point\n"%(no_pts,sur_idx))# put the point on surface
            else:
                if rz[i]>0:
                    fh.close()
                    raise ValueError("electrode z coordinate is greater than 0 in gmshWrap.py and you can't have refinement points above the surface!")
                fh.write("Point (%i) = {%.16f, %.16f, %.16f, %.16f};\n"%(no_pts, rx[i], ry[i], rz[i], rcl[i]))
                fh.write("Point{%i} In Volume{%i};//buried refinement point\n"%(no_pts,vol_idx))# put the point in volume 
        fh.write('//End mesh refinement points\n')
        
    if nfield>0: 
        nfield += 1 
        fh.write('//Refinement feild for fine mesh region\n')
        add_box_field(fh, nfield, min_x, max_x, min_y, max_y, 0, -fmd, 
                      cl*cl_factor, cln)
        fh.write('//Merge refinement feilds\n')
        set_fields(fh, nfield)
    fh.write('//End of script')
    
    fh.close()
    # print("writing .geo to file completed, save location:\n%s\n"%os.getcwd())
    return np.array(node_pos) 

#%% tank mesh (closed 3d box, no half-space)
def tank(elec=None, origin=None, dimension=[10.0,10.0,10.0],
         file_path='mesh3d.geo', cl=-1):
    """Created a .geo file for a closed 3D mesh (box or tank).

    Parameters
    ----------
    elec : list of list or 2D array, optional
        Electrodes positions given as X,Y,Z in each columns.
    origin : list of float, optional
        Origin of the corner where the mesh will be drawned from. If not
        provided and elec provided, the smaller elec position will be chosen.
    dimension : list of float, optional
        Dimension of the mesh in X,Y,Z from the corner origin.
        The default is [10,10,10].
    file_path : str, optional
        Path to the .geo file (with extension). The default is 'mesh3d.geo'.
    cl : float, optional
        Charateristic length of the electrode nodes. The default is -1.

    Returns
    -------
    list of float
        Indices of the electrode nodes.
        
    Notes:
    ------
    Code automatically identifies if electrodes lie on the tank's side

    """
    # checking arguments
    if file_path[-4:] != '.geo':
        file_path += '.geo'
    if elec is not None:
        if isinstance(elec, list):
            elec = np.array(elec)

    elec_x, elec_y, elec_z = elec[:,0], elec[:,1], elec[:,2]
    if origin is None:
        if elec is not None:
            origin = np.min(elec, axis=0)
        else:
            origin = [0,0,0]
            
    if cl < 0:
        if elec is None:
            cl = np.max(dimension)/10
        else:
            # compute distance matrix (useful for inferring cl and radius)
            dist = np.sort(np.unique(find_dist(elec_x,elec_y,elec_z)))[1:] # first is 0 ofc
            cl = np.min(dist)/2 # half of the minimum electrode spacing

    # compute corners
    min_x, max_x = np.sort([origin[0], origin[0]+dimension[0]])
    min_y, max_y = np.sort([origin[1], origin[1]+dimension[1]])
    min_z, max_z = np.sort([origin[2], origin[2]+dimension[2]])
    
    #check electrodes dont exceed the boundares of the tank
    if min(elec_x) < min_x or max(elec_x) > max_x:
        raise ValueError('X coordinate of an electrode exceeds tank bounds')
    if min(elec_y) < min_y or max(elec_y) > max_y:
        raise ValueError('Y coordinate of an electrode exceeds tank bounds')
    if min(elec_z) < min_z or max(elec_z) > max_z:
        raise ValueError('Z coordinate of an electrode exceeds tank bounds')
    
    def checkplane(p):
        #check if point lies on plane 
        S = np.cross(p[1] - p[0], p[2] - p[0])
        N = np.dot(S, p[0] - p[3])
        if N == 0:
            return True
        return False 
    
    # start writing the .geo
    with open(file_path, 'w') as fh:

        # headers
        fh.write("//3D tank mesh for ResIPy\n")
        fh.write("Mesh.Binary = 0;//specify we want ASCII format\n")
        fh.write("cl={:.2f}; //define characteristic length for electrode\n".format(cl))
        fh.write("cl2={:.2f}; //define characteristic length for tank mesh region\n".format(cl*1.2))
        
        fh.write("//Tank mesh region.\n")
        #add points for top of the box
        no_pts = 1
        nodes = []
        loop_pt_idx=[no_pts]
        nodes.append([max_x, max_y, max_z])
        fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl2};\n"%(no_pts, max_x, max_y, max_z))
        
        no_pts += 1
        loop_pt_idx.append(no_pts)
        nodes.append([max_x, min_y, max_z])
        fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl2};\n"%(no_pts, max_x, min_y, max_z))
        
        no_pts += 1
        loop_pt_idx.append(no_pts)
        nodes.append([min_x, min_y, max_z])
        fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl2};\n"%(no_pts, min_x, min_y, max_z))
        
        no_pts += 1
        loop_pt_idx.append(no_pts)
        nodes.append([min_x, max_y, max_z])
        fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl2};\n"%(no_pts, min_x, max_y, max_z))
        
        # add line loop
        no_lns = 0 
        for i in range(4):
            no_lns += 1 
            if i == 3:
                fh.write("Line(%i) = {%i,%i};\n"%(no_lns,loop_pt_idx[i],loop_pt_idx[0]))
            else:
                fh.write("Line(%i) = {%i,%i};\n"%(no_lns,loop_pt_idx[i],loop_pt_idx[i+1]))
        
        # bottom of the box  
        no_pts += 1
        loop2_pt_idx=[no_pts]
        nodes.append([max_x, max_y, min_z])
        fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl2};\n"%(no_pts, max_x, max_y, min_z))
        
        no_pts += 1
        loop2_pt_idx.append(no_pts)
        nodes.append([max_x, min_y, min_z])
        fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl2};\n"%(no_pts, max_x, min_y, min_z))
        
        no_pts += 1
        loop2_pt_idx.append(no_pts)
        nodes.append([min_x, min_y, min_z])
        fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl2};\n"%(no_pts, min_x, min_y, min_z))
        
        no_pts += 1
        loop2_pt_idx.append(no_pts)
        nodes.append([min_x, max_y, min_z])
        fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl2};\n"%(no_pts, min_x, max_y, min_z))
        
        # add line loops to connect the points just made 
        for i in range(4):
            no_lns += 1 
            if i == 3:
                fh.write("Line(%i) = {%i,%i};\n"%(no_lns,loop2_pt_idx[i],loop2_pt_idx[0]))
            else:
                fh.write("Line(%i) = {%i,%i};\n"%(no_lns,loop2_pt_idx[i],loop2_pt_idx[i+1]))
    
        # connect the top and bottom of the mesh 
        for i in range(4):
            no_lns += 1 
            fh.write("Line(%i) = {%i,%i};\n"%(no_lns,loop_pt_idx[i],loop2_pt_idx[i]))   
            
        fh.write("//End tank mesh region points\n" )
            
        
        # add relevant line loops for top surface      
        fh.write("Line Loop(1) = {1,2,3,4};\n") # fine mesh region top line loop         
        fh.write("Plane Surface(1) = {1};\n")
        
        # add sides and bottom surfaces of the box 
        fh.write("//Other tank mesh surfaces\n")
        fh.write("Line Loop(2) = {5,6,7,8};\n") # basal surface
        fh.write("Line Loop(3) = {5,-10,-1,9};\n") 
        fh.write("Line Loop(4) = {11,-6,-10,2};\n") 
        fh.write("Line Loop(5) = {11,7,-12,-3};\n") 
        fh.write("Line Loop(6) = {8,-9,-4,12};\n") 
        for i in range(2,7):
            fh.write("Plane Surface(%i) = {%i};\n"%(i,i))
        
        fh.write("Surface Loop (1) = {1,2,3,4,5,6};\n") # fine mesh region volume 
        fh.write("Volume(1) = {1}; //End tank mesh region surfaces.\n")
        
        surface_flag = np.zeros(elec.shape[0])
        surface_connection = [[1,2,3,4],
                              [5,6,7,8],
                              [6,5,2,1],
                              [6,7,3,2],
                              [7,8,3,4],
                              [8,5,4,1]]
        surface_connection = np.array(surface_connection) -1 
        node = np.array(nodes)
        
        for i in range(6):
            p = node[surface_connection[i,:]]
            for j in range(elec.shape[0]):
                p[-1,:] = elec[j,:]
                if checkplane(p):
                    surface_flag[j]=i+1
        
        node_pos = []
        if elec is not None:
            fh.write("//Electrode positions.\n")
            node_pos = [0]*len(elec_x)
            for i in range(len(elec_x)):
                no_pts += 1
                node_pos[i] = no_pts
                if surface_flag[i]>0:
                    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl};\n"%(no_pts, elec_x[i], elec_y[i], elec_z[i]))
                    fh.write("Point{%i} In Surface{%i}; //tank surface electrode\n"%(no_pts,surface_flag[i]))# put the point on surface
                else: # buried
                    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl};\n"%(no_pts, elec_x[i], elec_y[i], elec_z[i]))
                    fh.write("Point{%i} In Volume{1}; //buried (in tank) electrode\n"%(no_pts))# put the point in volume 
            
            fh.write("//End electrodes\n") 

    return node_pos


#%% general column mesh (built with prisms)
def prism(electrodes, poly=None, z_lim=None, radius=None,
          file_path='prism_mesh.geo', cl=-1, elemz=4):
    """Make a prism mesh.
    
    Parameters
    ------------
    electrodes: list of array likes
        first column/list is the x coordinates of electrode positions and so on ... 
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
    if file_path.find('.geo')==-1:
        file_path=file_path+'.geo'#add file extension if not specified already
    
    elec_x = electrodes[0]
    elec_y = electrodes[1]
    elec_z = electrodes[2]
    #uni_z = np.unique(elec_z) # unique z values 
    
    if z_lim is None:
#        print('zlim not passed')
        z_lim = [min(elec_z),max(elec_z)]

    if radius is None:
        radius  = max(electrodes[0])
    if cl == -1:
        dist_sort = np.unique(find_dist(elec_x,elec_y,elec_z))
        cl = dist_sort[1]/2 # characteristic length is 1/2 the minimum electrode distance   
    
    num_elec = len(elec_z)

    fh = open(file_path,'w')
    fh.write("// ResIPy column (or prism) mesh script\n")
    fh.write('SetFactory("OpenCASCADE");\n')
    fh.write("Mesh.Binary = 0;//specify we want ASCII format\n")    
    fh.write("cl=%f;\n"%cl)
    
    x = []
    y = [] # ignore repeated vertices 
    
    if poly is None:
        fh.write("// Make a circle which is extruded\n")
        fh.write("Point (1) = {0,0,%f,cl};\n"%(z_lim[0]))
        x.append(0);y.append(0)
        fh.write("Point (2) = {%f,0,%f,cl};\n"%(radius,z_lim[0]))
        x.append(radius);y.append(0)
        fh.write("Point (3) = {0,%f,%f,cl};\n"%(radius,z_lim[0]))
        x.append(0);y.append(radius)
        fh.write("Point (4) = {0,-%f,%f,cl};\n"%(radius,z_lim[0]))
        x.append(0);y.append(-radius)
        fh.write("Point (5) = {-%f,0,%f,cl};\n"%(radius,z_lim[0]))
        x.append(-radius);y.append(0)
        fh.write("Circle (1) = {2,1,3};\n")
        fh.write("Circle (2) = {3,1,4};\n")
        fh.write("Circle (3) = {4,1,5};\n")
        fh.write("Circle (4) = {5,1,2};\n")

        fh.write("Line Loop(1) = {3, 4, 1, 2};\n")
        fh.write("Plane Surface(1) = {1};\n")
        fh.write("Point {%i} In Surface {1};\n\n"%1)
        edges = 3
        pt_no = 6
        
        #work out which electrodes are inside or on the end of the column 
        dist = np.sqrt(np.array(elec_x)**2 + np.array(elec_y)**2)
        inside = dist < radius 
    else:
        poly_x = poly[0]
        poly_y = poly[1]
        pt_no = 0
        #make a node for each vertex in the polygon 
        for i in range(len(poly_x)):
            pt_no += 1
            fh.write("Point (%i) = {%f,%f,%f,cl};//polygon point\n"%(pt_no,poly_x[i],poly_y[i],z_lim[0]))
        # connect the nodes with lines 
        ln_no = 0
        for i in range(pt_no-1):
            ln_no += 1
            fh.write("Line (%i) = {%i,%i};//polygon line\n"%(ln_no,i+1,i+2))
        ln_no+=1
        fh.write("Line (%i) = {%i,%i};//closing line\n"%(ln_no,pt_no,1))
        edges = ln_no
        #form a line loop and surface
        fh.write("Line Loop(1) = {")
        [fh.write('%i,'%(i+1)) for i in range(ln_no-1)]
        fh.write("%i};\n"%ln_no)
        fh.write("Plane Surface(1) = {1};\n\n")
        
        #work out which electrodes are inside or on the end of the column 
        path = mpath.Path(np.array([poly[0],poly[1]]).T)
        inside = path.contains_points(np.array([elec_x,elec_y]).T)
        
    
    if all(inside) == False:
        fh.write("//Points inside the column\n")
        x = []
        y = [] # ignore repeated vertices 
        for i in range(num_elec):
            if inside[i] and elec_x[i] not in x and elec_y[i] not in y:
                pt_no+=1
                fh.write("Point (%i) = {%f,%f,%f,cl};\n"%(pt_no,elec_x[i],elec_y[i],z_lim[0]))
                fh.write("Point {%i} In Surface {1};\n"%pt_no)
                x.append(elec_x[i])
                y.append(elec_y[i])
            
    surface = 1
    
    if min(elec_z)<z_lim[0] or max(elec_z) > z_lim[1]:
        fh.close()
        raise ValueError ("Z coordinate of column electrodes is out of bounds")
    
    allz = np.append(elec_z,z_lim)
    
    uni_z = np.unique(allz)
    #extrude surfaces 
    fh.write("//Extrude planes in between each electrode.\n") 
    seg = 0 # segment number 
        
    for i in range(0,len(uni_z)-1):
        #work out the amount to extrude 
        diff = uni_z[i+1] - uni_z[i]    
        seg += 1          
        fh.write('seg%i = Extrude {0, 0, %f} {Surface{%i}; Layers{%i}; Recombine;};'%(seg,diff,surface,elemz))
        surface += edges + 1
        if i==0:
            fh.write('//Base of column\n')
        elif i == len(uni_z)-2:
            fh.write('//top of column\n')
        else:
            fh.write('\n')
        
    fh.write('Physical Volume(1) = {')
    for i in range(seg-1):
        fh.write('seg%i[1],'%(i+1))
    fh.write('seg%i[1]};'%(seg))
        
    fh.close()
    
#%% cylinder mesh 
def cylinder(electrodes, zlim=None, radius=None,
             file_path='cylinder_mesh.geo', cl=-1, elemz=4,cl_factor=2, 
             add_refine=True):
    """Make a cylinderical tetrahedral mesh.
    
    Parameters
    ------------
    electrodes: list of array likes
        first column/list is the x coordinates of electrode positions and so on ... 
        geometry assumes centre of column base is at the origin. 
    poly: list, tuple, optional 
        Describes polygon where the argument is 2 by 1 tuple/list. Each entry is the polygon 
        x and y coordinates, ie (poly_x, poly_y)
    zlim: list, tuple, optional 
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
        Inactive parameter, left for compatiblity with prism mesh function. 
    cl_factor: float, optional 
        Factor by which to grow elements away from electrodes on the column. 
        Default is 2. 
    add_refine: bool, optional
        Add in points to refine mesh around electrodes only, default is True. 
    """
    def cdist(x0,y0,x1,y1):
        return np.sqrt((x0-x1)**2 + (y0-y1)**2)
    # find unique x y 
    if file_path.find('.geo')==-1:
        file_path=file_path+'.geo'#add file extension if not specified already
    pt_no = 1 
    elec_x = np.array(electrodes[0])
    elec_y = np.array(electrodes[1])
    elec_z = np.array(electrodes[2])
    nelec = len(elec_x)
    
    #work out which electrodes are inside or on the end of the column 
    dist = np.sqrt(np.array(elec_x)**2 + np.array(elec_y)**2)
    # we assume the origin is at 0,0 
    if radius is None:
        radius = np.max(dist) 

    diff = np.abs(dist - radius)
    inside = diff > 1e8 # electrodes on inside of column radius 
    
    rx = np.linspace(-radius,radius,1000000)
    ry = np.sqrt(radius**2 - rx**2)
    cx = np.append(rx,rx)
    cy = np.append(ry,-ry)
    
    # do correction for electrodes not quite exactly on radius
    for i in range(nelec):
        if not inside[i] and dist[i] != radius:
            dx = elec_x[i] - cx 
            dy = elec_y[i] - cy 
            dd = np.sqrt((dx**2)+(dy**2))
            di = np.argmin(dd)
            elec_x[i] = cx[di]
            elec_y[i] = cy[di]
            print('Electrode %i not quite right, adjusting...'%i)
    
    # find unique xy positions
    uni_x = []
    uni_y = [] 
    uni_flag = []
    uni_ref = []
    curc_x = [] # circumfrence points 
    curc_y = []

    for i in range(nelec):
        if len(uni_x) == 0: # first point is here 
            uni_x.append(elec_x[i])
            uni_y.append(elec_y[i])
            uni_flag.append(inside[i])
            uni_ref.append(False)
            if not inside[i]:
                curc_x.append(elec_x[i])
                curc_y.append(elec_y[i])
        dist = [cdist(uni_x[j],uni_y[j],elec_x[i],elec_y[i]) for j in range(len(uni_x))]
        if min(dist) > 0: # if distance greater than 0 then the point has not been seen before 
            uni_x.append(elec_x[i])
            uni_y.append(elec_y[i])
            uni_flag.append(inside[i]) # flag for if unique point on outside or inside of column 
            uni_ref.append(False) # flag if refinement point 
            if not inside[i]:
                curc_x.append(elec_x[i])
                curc_y.append(elec_y[i])
            
    # check for at least 4 unique points on column circumfrence on all quadrants of the column 
    quad_check = False 
    if len(curc_x) < 4: 
        quad_check = True 
    
    # now check all quadrants are occupied, if not add extra points in all quads 
    print(curc_x, curc_y)
    if min(curc_x) <= 0 and max(curc_x) <= 0: 
        quad_check = True 
    if min(curc_x) >= 0 and max(curc_x) >= 0: 
        quad_check = True 
    if min(curc_y) <= 0 and max(curc_y) <= 0: 
        quad_check = True 
    if min(curc_y) >= 0 and max(curc_y) >= 0: 
        quad_check = True 
        
    if quad_check: 
        #make constructor points for making a basic column 
        cons_x = [0,radius,0,-radius]
        cons_y = [radius,0,-radius,0]
        for i in range(4):
            dist = [cdist(uni_x[j],uni_y[j],cons_x[i],cons_y[i]) for j in range(len(uni_x))]
            if min(dist) > 0: # if distance greater than 0 then the point has not been seen before 
                uni_x.append(cons_x[i])
                uni_y.append(cons_y[i])
                uni_flag.append(False) # flag for if unique point on outside or inside of column 
                uni_ref.append(True) # flag if refinement point 
                curc_x.append(cons_x[i])
                curc_y.append(cons_y[i])
        
    bears = [bearing(curc_x[i],curc_y[i]) for i in range(len(curc_x))]
    bear_sort = np.argsort(bears)
    curc_x = np.array(curc_x)[bear_sort]
    curc_y = np.array(curc_y)[bear_sort]
    
    if add_refine:        
        # were gonna walk around the circle and pick out coordinates for refinement points 
        for i in range(len(curc_x)):
            if i == len(curc_x)-1:
                dx = curc_x[0] - curc_x[i]
                dy = curc_y[0] - curc_y[i]
            else:
                dx = curc_x[i+1] - curc_x[i]
                dy = curc_y[i+1] - curc_y[i]
            d = np.sqrt(dx**2 + dy**2)
            if d < 0.5*cl: 
                # skip if inter electrode distance smaller than 50% characteristic length 
                continue 
            # mid point of circumfrence points 
            mx = curc_x[i] + (dx/2)
            my = curc_y[i] + (dy/2)

            dist = np.sqrt((cx-mx)**2 + (cy-my)**2)
            mi = np.argmin(dist)
            uni_x.append(cx[mi])
            uni_y.append(cy[mi])
            uni_flag.append(False)
            uni_ref.append(True)
            
    # find bearings for unique xy points 
    bears = [bearing(uni_x[i],uni_y[i]) for i in range(len(uni_x))]
    bear_sort = np.argsort(bears)
    # sort by bearing 
    uni_x = np.array(uni_x)[bear_sort].tolist()
    uni_y = np.array(uni_y)[bear_sort].tolist()
    uni_flag = np.array(uni_flag,dtype=bool)[bear_sort].tolist()
    uni_ref = np.array(uni_ref,dtype=bool)[bear_sort]
            
    if zlim is None:
        zlim = [min(elec_z),max(elec_z)]
        
    if cl == -1:
        dist_sort = np.unique(find_dist(elec_x,elec_y,elec_z))
        cl = dist_sort[1]/2 # characteristic length is 1/2 the minimum electrode distance   
    
    num_elec = len(elec_z)
    
    fh = open(file_path,'w')
    fh.write("// ResIPy tetrahedral column mesh script\n")
    fh.write('SetFactory("OpenCASCADE");\n')
    fh.write("Mesh.Binary = 0;//specify we want ASCII format\n")    
    fh.write("cl=%f;\n"%cl)
    fh.write("cl_fac=%f;\n"%cl_factor)
    
    fh.write("// Make a circle at lowest z limit which is extruded\n")
    if any(inside):
        fh.write("Point (%i) = {0,0,%f,cl};\n"%(pt_no,zlim[0]))
    else:
        fh.write("Point (%i) = {0,0,%f,cl*cl_fac};\n"%(pt_no,zlim[0]))

    ## go around clockwise  
    for i in range(len(uni_x)):
        if not uni_flag[i]: 
            pt_no += 1 
            fh.write("Point (%i) = {%32.16f,%32.16f,%32.16f,cl"%(pt_no,uni_x[i],uni_y[i],zlim[0]))
            if uni_ref[i]:
                fh.write("*cl_fac};//column circumfrence point(coarse)\n")
            else:
                fh.write("};//column circumfrence point (fine)\n")
    
    ncirc_points = pt_no - 1     
        
    circle_no = 0
    for i in range(1,pt_no):
        circle_no += 1 
        if i+1 == pt_no:
            fh.write("Circle (%i) = {%i,1,%i};//circumfrence segment(last one)\n"%(circle_no,i+1,2))
        else:
            fh.write("Circle (%i) = {%i,1,%i};//circumfrence segment\n"%(circle_no,i+1,i+2))
    # edges = circle_no 
         
    fh.write("Line Loop(1) = {1")
    for i in range(1,circle_no):
        fh.write(",%i"%(i+1))
    fh.write("};//Base of column\n")
    
    # make a circle in gmsh 
    fh.write("Plane Surface(1) = {1};\n")
    fh.write("Point {%i} In Surface {1};\n\n"%1)

    fh.write("//Points on surface of the column\n")
    for i in range(len(uni_x)):
        if uni_x[i] == 0 and uni_y[i]==0:
            continue # continue as this point is at origin 
        if uni_flag[i]: # check if point inside of column radius, hence add it to column 
            pt_no+=1
            fh.write("Point (%i) = {%32.16f,%32.16f,%32.16f,cl};\n"%(pt_no,uni_x[i],uni_y[i],zlim[0]))
            fh.write("Point {%i} In Surface {1};\n"%pt_no)

    fh.write("//End points on surface of column\n")
    surface = 1
    
    if min(elec_z)<zlim[0] or max(elec_z) > zlim[1]:
        fh.close()
        raise ValueError ("Z coordinate of column electrodes is out of bounds")
    
    allz = np.append(elec_z,zlim)
    
    uni_z = np.unique(allz)
    max_diff = np.max(np.abs(np.diff(uni_z)))
    refine_z = [] 
    if add_refine: # triggers if refinement requested and no internal electrodes 
        for i in range(1,len(uni_z)):
            diff = abs(uni_z[i] - uni_z[i-1]) 
            if diff<(0.2*max_diff):
                continue 
            refine_z.append(uni_z[i-1] + (diff/2))

    uni_z = np.unique(np.append(uni_z,refine_z))
    
    #extrude surfaces 
    fh.write("//Extrude planes in between each electrode.\n") 
    seg = 0 # segment number 
    nseg = len(uni_z)-1 # number of segments 
    point_cache = np.zeros((nseg,ncirc_points),dtype=int)
    refinement = np.zeros(nseg,dtype=int)
    
    for i in range(nseg):
        #work out the amount to extrude 
        diff = uni_z[i+1] - uni_z[i]    
        seg += 1          
        fh.write('seg%i = Extrude {0, 0, %f} {Surface{%i};Recombine;};'%(seg,diff,surface))
        surface += ncirc_points + 1 
        if uni_z[i+1] in elec_z:
            fh.write('//Electrode plane')
        else:
            fh.write('//Refinement plane')
            refinement[i] = 1 
        if i == nseg-1:
            fh.write('(Also Top of column)\n')
        else:
            fh.write('\n')
        point_cache[i,:] = [pt_no + 1 + j for j in range(ncirc_points)] 
        pt_no += ncirc_points # add points in circle for each extruded segment 

        
    # delete extruded volumes 
    fh.write('Delete {')
    for i in range(seg):
        fh.write('Volume{%i};'%(i+1))
    fh.write('}\n')

    top_surf = 1 # number of circle surface at top of colum 
    base_surf = surface # number of circle surface at base of column 

    c = 0 
    fh.write('Surface Loop(%i) = {%i,'%(seg+1,top_surf))
    for i in range(1,base_surf):
        if c == ncirc_points:
            c = 0 
            continue 
        fh.write('%i,'%(i+1))
        c+=1 
    fh.write('%i};\n'%base_surf)
    fh.write('Volume(1) = {%i};//make surface loop a volume\n'%(seg+1)) 
    
    # check for electrodes in the column 

    fh.write("//Points inside the column\n")
    for i in range(num_elec):
        if elec_z[i] == zlim[0] or elec_z[i] == zlim[1]:
            continue # should be included in the surface already 
        if inside[i]:
            pt_no+=1
            fh.write("Point (%i) = {%f,%f,%f,cl};//Internal electrode\n"%(pt_no,elec_x[i],elec_y[i],elec_z[i]))
            fh.write("Point {%i} In Volume {1};\n"%pt_no)
            
    for i in range(1,len(uni_z)):
        if all(inside==False):
            continue # skip if internal electrodes present 
        z = uni_z[i]
        diff = abs(uni_z[i] - uni_z[i-1])
        if z== zlim[0] or z == zlim[1]:
            continue # should be included in the surface already 
        if diff<(0.2*max_diff):
            continue 
        pt_no +=1 
        fh.write("Point (%i) = {%f,%f,%f,cl*cl_fac};//Internal refinement\n"%(pt_no,0,0,z))
        fh.write("Point {%i} In Volume {1};\n"%pt_no)
        
    fh.write("//End Points inside the column\n")
    fh.write('Physical Volume(1) = {1};//make column one physical volume\n')
    
    # make base of column points bigger if not actually in line with electrodes 
    if zlim[0] != min(elec_z):
        points = np.arange(ncirc_points)+2 
        point_cache = np.vstack([points,point_cache])
        refinement = np.append(np.ones(1,dtype=int),refinement)
        fh.write("Characteristic Length {1} = cl*cl_fac;\n")
    
    # grow mesh elements away from electrodes
    fh.write("//Grow mesh elements away from electrodes\n")
    for i in range(len(refinement)):
        if refinement[i] == 1:
            points = point_cache[i,:]
            fh.write("Characteristic Length {")
            for j in range(ncirc_points):
                if j > 0:
                    fh.write(',')
                fh.write('%i'%points[j])
            fh.write("} = cl*cl_fac;\n")
    # add a final centre point 
    # surface += 1 
    # pt_no += 1
    # fh.write('//Top of column middle point\n')
    # fh.write("Point (%i) = {%f,%f,%f,cl};\n"%(pt_no,0,0,zlim[1]))
    # fh.write("Point {%i} In Surface {%i};\n"%(pt_no,surface))
    
    fh.write('//End of script\n')
    fh.close()


#%% Cylinder mesh using tetrahedra
def cylinder_mesh_old(electrodes, zlim=None, radius=None, 
                  file_path='cylinder_mesh.geo', cl=-1, finer=4):
    """Make a cylinder mesh using tetrahedra.
    
    Parameters
    ------------
    electrodes : list of array_like or nx3 array
        First column/list is the x coordinates of electrode positions and so on ... 
    zlim : list, tuple, optional 
        Bottom and top z coordinate of column, in the form (min(z),max(z))
    radius: float, optional 
        Radius of column.
    file_path : string, optional 
        Name of the generated gmsh file (can include file path also).
    cl : float, optional
        Characteristic length, essentially describes how big the nodes 
        associated elements will be. Usually no bigger than 5. If set as -1 (default)
        a characteristic length 1/4 the minimum electrode spacing is computed. 
    finer : int, optional
        Number of line between two consecutive electrodes to approximate the
        circle shape.
    """
    # check the type of electrodes argument
    if isinstance(electrodes, list):
        elec = np.array(electrodes).T
    else:
        elec = electrodes
        
    # compute number of electrode per ring
    nepr = np.vstack(list({tuple(e) for e in elec[:, :2]})).shape[0] 
    
    # compute distance matrix (useful for inferring cl and radius)
    def cdist(a):
        z = np.array([complex(x[0], x[1]) for x in a])
        return np.abs(z[...,np.newaxis]-z)

    uz = np.unique(elec[:,2])
    ie = elec[:,2] == uz[0]
    dist = cdist(elec[ie,:2]).flatten()
    
    content = ''
    
    # characteristic length
    if cl == -1:
        cl = np.min(dist[dist > 0])/3
    content = content + 'cl = {:f};\n'.format(cl)
    
    c = 1 # entity count
    # add electrodes
    content += '// Electrodes locations\n'
    cpts = []
    for i, pts in enumerate(elec):
        content += 'Point({:d}) = {{{:f}, {:f}, {:f}, cl}};\n'.format(c, *pts)
        cpts.append(c)
        c += 1
    
    # add file extension if not specified already
    if file_path[-4:] != '.geo':
        file_path = file_path + '.geo'
        
    # compute radius of not provided
    if radius is None:
        radius = np.max(dist)/2
        
    # compute default zlim if not provided
    if zlim is None:
        zlim = [np.max(elec[:,2]), np.min(elec[:,2])]
        
    # compute angles and stages
    angles = np.linspace(0, 2*np.pi, nepr*finer+1)[:-1]
    celec = np.c_[radius*np.cos(angles), radius*np.sin(angles)]
    stages = [np.max(zlim), np.min(zlim)]
    ringPerStage = [len(np.unique(elec[:,2]))]
    
    # add additional points to top and bottom ring to make it more circular
    cstg = []
    for stage in stages:
        tmp = []
        for i, a in enumerate(celec):
            content += 'Point({:d}) = {{{:f}, {:f}, {:f}, cl}};\n'.format(c, a[0], a[1], stage)
            tmp.append(c)
            c += 1
        tmp = tmp + [tmp[0]]
        cstg.append(tmp)
    
    # horizontal lines
    lhor = []
    for a in cstg:
        tmp = []
        for i in range(nepr*finer):
              content += 'Line({:d}) = {{{:d}, {:d}}};\n'.format(c, a[i], a[i+1])
              tmp.append(c)
              c += 1
        tmp = tmp + [tmp[0]]
        lhor.append(tmp)
    
    # vertical lines
    lverAll = []
    for l, ring in enumerate(ringPerStage):
        lver = []
        atop = cstg[l]
        abot = cstg[l+1]
        for i, ce in enumerate(celec):
            # is there any electrode aligned with this node?
            ie = np.isclose(elec[:, 0], ce[0]) & np.isclose(elec[:, 1], ce[1])
            
            if np.sum(ie) > 0: # electrode line, let's connect
                # b = int(i/finer + 1)
                # offset = 0 if l == 0 else np.sum(np.array(ringPerStage)[:l])*nepr
                # print(b, offset, i, finer, atop, abot, lver, ring)
                # el = np.arange(offset + b, offset + b + ring*nepr, nepr).astype(int).tolist()
                iel = np.where(ie)[0] + 1  # first elec is 1, not 0 in .geo file
                isort = np.argsort(elec[ie, 2])[::-1]  # from high z to low z
                el = list(iel[isort])
                a = [atop[i]] + el + [abot[i]]  # list concatenation
                tmp = []
                for j in range(len(a)-1):
                    content += 'Line({:d}) = {{{:d}, {:d}}};\n'.format(c, a[j], a[j+1])
                    tmp.append(c)
                    c += 1
                lver.append(tmp)
            else:
                content += 'Line({:d}) = {{{:d}, {:d}}};\n'.format(c, atop[i], abot[i])
                lver.append([c])
                c += 1
        lver = lver + [lver[0]]
        lverAll.append(lver)
    
    # line loop
    sver = []
    for l, lver in enumerate(lverAll):
        ll = []
        ltop = lhor[l]
        lbot = lhor[l+1]
        for i in range(nepr*finer):
            a = str(lver[i][0]) if len(lver[i]) == 1 else ' ,'.join(np.array(lver[i]).astype(int).astype(str))
            b = str(-lver[i+1][0]) if len(lver[i+1]) == 1 else ' ,'.join((-np.array(lver[i+1]))[::-1].astype(int).astype(str))
            content += 'Line Loop({:d}) = {{{:d}, {:s}, {:d}, {:s}}};\n'.format(
                c, -ltop[i], a, lbot[i], b)
            ll.append(c)
            c += 1
        sver.append(ll)
    
    # horizontal line loop
    shor = []
    for l in lhor:
        # n = nepr*finer
        content += 'Line Loop({:d}) = {{'.format(c) + ' ,'.join(np.array(l)[:-1].astype(int).astype(str)) + '};\n'
        shor.append(c)
        c += 1
    
    # surface
    slv = []
    for i, ll in enumerate(sver):
        tmp = []
        for l in ll:
            content += 'Plane Surface({:d}) = {{{:d}}};\n'.format(c, l)
            tmp.append(c)
            c += 1
        slv.append(tmp)
    
    slh = []
    for i, l in enumerate(shor):
        content += 'Plane Surface({:d}) = {{{:d}}};\n'.format(c, l)
        slh.append(c)
        c += 1
        
    # surface loop and volume
    for i in range(len(stages)-1):
        a = [slh[i]] + slv[i] + [slh[i+1]]
        content += 'Surface Loop({:d}) = {{'.format(c) + ' ,'.join(np.array(a).astype(str)) + '};\n'
        content += 'Volume({:d}) = {{{:d}}};'.format(i, c)
        c += 1
        
    # write content to file
    with open(file_path, 'w') as f:
        f.write("Mesh.Binary = 0;//specify we want ASCII format\n")
        f.write(content)


# test
# radius = 6.5/2 # cm
# angles = np.linspace(0, 2*np.pi, 13)[:-1][::6] # radian
# celec = np.c_[radius*np.cos(angles), radius*np.sin(angles)]
# elec2 = np.c_[np.tile(celec.T, 8).T, np.repeat(6.5+np.arange(0, 8*5.55, 5.55)[::-1], len(angles))]
# cylinder_mesh(elec, file_path='invdir/mesh_cylinder.geo',
#               zlim=[0, 48])


    
#%% Make an arbitary 2D shape 
def mesh2d(electrodes, poly=None, origin=(0,0), radius=None,
           file_path='2d_mesh.geo', cl=-1):
    """Make generic 2d mesh.
    
    Parameters
    ------------
    electrodes: list of array likes
        first column/list is the x coordinates of electrode positions and so on ... 
    poly: list, tuple, optional 
        Describes polygon where the argument is 2 by 1 tuple/list. Each entry is the polygon 
        x and y coordinates, ie (poly_x, poly_y)
    origin: list, tuple, optional 
        2 by 1 array of the x and y coordinates of the circle origin, only used if poly is None. 
    radius: float, optional 
        radius of column
    file_path: string, optional 
        name of the generated gmsh file (can include file path also) (optional)
    cl: float, optional
        characteristic length (optional), essentially describes how big the nodes 
        assocaited elements will be. Usually no bigger than 5. If set as -1 (default)
        a characteristic length 1/4 the minimum electrode spacing is computed.
    """
    if file_path.find('.geo')==-1:
        file_path=file_path+'.geo'#add file extension if not specified already
    
    elec_x = np.array(electrodes[0])
    elec_y = np.array(electrodes[1])
    elec_z = np.array(len(elec_x)*[0])
    
    ox = origin[0]
    oy = origin[1]

    if radius is None:
        radius  = max(electrodes[0])
    if cl == -1:
        dist_sort = np.unique(find_dist(elec_x,elec_y,elec_z))
        cl = dist_sort[1]/2 # characteristic length is 1/2 the minimum electrode distance   
        
    r = radius
    
    num_elec = len(elec_z)

    fh = open(file_path,'w')
    fh.write("// ResIPy 2d shape mesh script\n")
    fh.write("Mesh.Binary = 0;//specify we want ASCII format\n")
    fh.write('SetFactory("OpenCASCADE");\n')
    fh.write("cl=%f;\n"%cl)
    
    if poly is None:
        fh.write("// Make a circle\n")
        fh.write("Point (1) = {%f,%f,%f,cl};\n"%(ox,oy,0))
        fh.write("Point (2) = {%f,%f,%f,cl};\n"%(ox+r,oy,0))
        fh.write("Point (3) = {%f,%f,%f,cl};\n"%(ox,oy-r,0))
        fh.write("Point (4) = {%f,%f,%f,cl};\n"%(ox-r,oy,0))
        fh.write("Point (5) = {%f,%f,%f,cl};\n"%(ox,oy+r,0))
        fh.write("Circle (1) = {2,1,3};\n")
        fh.write("Circle (2) = {3,1,4};\n")
        fh.write("Circle (3) = {4,1,5};\n")
        fh.write("Circle (4) = {5,1,2};\n")

        fh.write("Line Loop(1) = {3, 4, 1, 2};\n")
        fh.write("Plane Surface(1) = {1};\n")
        fh.write("Point {1} In Surface {1};\n\n")
        pt_no = 6
        
        #work out which electrodes are inside or on the end of the column 
        dist = np.sqrt((elec_x-ox)**2 + (elec_y-oy)**2)
        inside = dist+0.00001 < radius 
        x = [ox,ox-r,ox+r]
        y = [oy,oy-r,oy+r]
    else:
        poly_x = poly[0]
        poly_y = poly[1]
        pt_no = 0
        #make a node for each vertex in the polygon 
        for i in range(len(poly_x)):
            pt_no += 1
            fh.write("Point (%i) = {%f,%f,%f,cl};//polygon point\n"%(pt_no,poly_x[i],poly_y[i],0))
        # connect the nodes with lines 
        ln_no = 0
        for i in range(pt_no-1):
            ln_no += 1
            fh.write("Line (%i) = {%i,%i};//polygon line\n"%(ln_no,i+1,i+2))
        ln_no+=1
        fh.write("Line (%i) = {%i,%i};//closing line\n"%(ln_no,pt_no,1))
        #form a line loop and surface
        fh.write("Line Loop(1) = {")
        [fh.write('%i,'%(i+1)) for i in range(ln_no-1)]
        fh.write("%i};\n"%ln_no)
        fh.write("Plane Surface(1) = {1};\n\n")
        
        #work out which electrodes are inside or on the end of the column 
        path = mpath.Path(np.array(poly[0],poly[1]).T)
        inside = path.contains_points(np.array([elec_x,elec_y]).T)
        x = poly_x
        y = poly_y
    
    if all(inside) == False:
        fh.write("//Points inside the shape\n")
         # ignore repeated vertices 
        for i in range(num_elec):
            if inside[i] and elec_x[i] not in x or elec_y[i] not in y:
                pt_no+=1
                fh.write("Point (%i) = {%f,%f,%f,cl};\n"%(pt_no,elec_x[i],elec_y[i],0))
                fh.write("Point {%i} In Surface {1};\n"%pt_no)
                x.append(elec_x[i])
                y.append(elec_y[i])
    
    fh.write("//End\n")    
    fh.close()
