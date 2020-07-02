Desc = """ 
Convert any .msh file into a .dat file for R2 and R3t. Code is adapted from 
the ResIPy package, an all in one python wrapper for R2, cR2, R3t and cR3t. 
The package even includes a graphical user interface. It can be found at 
https://gitlab.com/hkex/pyr2

Author: Jimmy Boyd - jamyd91@bgs.ac.uk
Date: 3rd of December 2019
"""
note = """
Notes: 
    -The code is designed to work with gmsh version 3.06 and version 4+ 

    -Code will accept command line options for the filepaths, however they are 
    not needed as file dailogue boxes will be opened to choose/save the relevant
    files.
    
    -If you have made "zones" in the mesh using physical entities, then pass 
    "--zones True", otherwise they will be ignored. Often it helps to zone the 
    mesh for construction purposes within gmsh, but not for inversion.

"""
#import modules 
import numpy as np 
import warnings, time, argparse
import tkinter as tk
from tkinter import filedialog

#%% Define relevant functions
def ui_open_file(title = "Please select gmsh mesh file",file_type=["mesh file","msh"]):
    root=tk.Tk()
    root.withdraw()
    file_path=filedialog.askopenfilename(title=title,filetypes=((file_type[0],"*."+file_type[1]),("all files","*.*")))
    return file_path

def ui_new_file(title = "Please save mesh file as",file_type=["R2/R3t mesh file","dat"]):
    root=tk.Tk()
    root.withdraw()
    file_path=filedialog.asksaveasfilename(title=title,filetypes=((file_type[0],"*."+file_type[1]),("all files","*.*")))
    return file_path

def ccw(p,q,r):#code expects points as p=(x,y) and so on ... 
    """Checks if points in a triangle are ordered counter clockwise. When using R2,
    mesh nodes should be given in a counter clockwise order otherwise you'll get negative 
    apparent resistivities. 
    
    Parameters
    ----------
    p - tuple, list,
        The x y coordinates of the point/vertex ie. (x,y) in desired order
    q - " 
    r - " 
    
    Returns
    ----------
    0 if colinear points, 1 if counter clockwise order, 2 if points are ordered clockwise
    """
    val=((q[1]-p[1])*(r[0]-q[0]))-((q[0]-p[0])*(r[1]-q[1]))
    if val==0:
        return 0 # lines are colinear
    elif val>0:
        return 1 # points are oreintated counter clockwise 
    elif val<0:
        return 2 # points are counter clockwise
    
def check_tetra(x,y,z):
    """Check if 4 points in a tetrahedra are ordered counter clocwise
    Parameters
    -----------
    x : array like 
        x coordinates of tetrahedral cell, 4 by 1 array. 
    y : array like 
        y coordinates of tetrahedral cell, 4 by 1 array. 
    z : array like 
        z coordinates of tetrahedral cell, 4 by 1 array.
    
    Returns
    ----------
    0 if coplanar points, 1 if clockwise order, 2 if points are counter clockwise
    """

    p = np.array([x,y,z]).T
    
    S = np.cross(p[1] - p[0], p[2] - p[0])
    N = np.dot(S, p[0] - p[3])
    
    if N>0:
        return 1 # points are clockwise
    elif N<0:
        return 2 # points are counter clockwise
    else:
        return 0 # something dogdey has happened as all points are on the same plane 

# parse a .msh file
def msh_parse(file_path):
    """Converts a gmsh mesh file into 
    
    Parameters
    ----------
    file_path: string
        file path to mesh file. note that a error will occur if the file format is not as expected
   
    Returns
    ----------
    Mesh class
    """
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
    
    print('Reading %s'%file_path)
    
    if mesh_format == '2.2 0 8':
        print('Gmsh version == 3.x')
        gmshV = 3 # assume its gmsh version 3.06
    else:
        print('Gmsh version == 4.x')
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
    print('reading node coordinates...')
    #read in number of nodes - at line 5
    #allocate lists for node numbers and coordinates
    node_num=[0]*no_nodes
    nodex=[0]*no_nodes
    nodey=[0]*no_nodes
    nodez=[0]*no_nodes
    #read in node information
    node_idx = 1
    for i in range(node_start, node_end):
        line_info=dump[i].split()
        #convert string info into floats
        line=[float(k) for k in line_info]
        previous_node = node_idx 
        node_idx = int(line[0])
        if not node_idx < previous_node: # else then its some weird flag used by gmsh, ignore
            node_num[node_idx-1]=node_idx
            nodex[node_idx-1]=line[1]
            nodey[node_idx-1]=line[2]
            nodez[node_idx-1]=line[3]
    
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
    print('Determining element type...',end='') # this depends a bit on the version of gmsh 
    for i in range(element_start,element_end):
        line = dump[i].split()
        if gmshV == 3:
            elm_type.append(int(line[1]))
        else:
            elm_type.append(len(line) - 1)
    
    #if using gmsh version 4.x filter out the flag lines 
    if gmshV == 4:
        flagline = []
        flagidx = []
        flaglines = int(dump[element_start-1].split()[0])
    
        line = dump[element_start].split()
        flagline.append(dump[element_start])
        flagidx.append(element_start)
        skip = int(line[-1])
        for i in range(flaglines-1): 
            line = dump[element_start+skip+1].split()
            flagline.append(dump[element_start+skip+1])
            flagidx.append(element_start+skip+1)
            skip = int(line[-1])
    
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
        print('Triangle')
        con_matrix = [[],[],[]]
    elif npere == 4: 
        print('Tetrahedra')
        con_matrix = [[],[],[],[]]
    elif npere == 6:
        print('Prism')
        con_matrix = [[],[],[],[],[],[]]
    else:
        raise ValueError('Cannot parse mesh becuase the relevant cell types cannot be found')
        
    print('Reading connection matrix...')
        
    phys = 0
    for i in range(element_start,element_end):
        splitted = dump[i].split()
        line=[int(k) for k in splitted]
        if gmshV == 3: 
            elmV = line[1]
            if npere == 3:
                elmV+=1 # in this format the flag for triangle elements is 2 so add 1
            st = 5
        else:
            elmV = len(line) - 1
            st = 1
            if i in flagidx:
                phys += 1 # physical entity incriments with flag lines 
                elmV = -1 # then ignore this line as it's a flag
        
        #convert string info into floats and cache data
        if npere == elmV:
            nat_elm_num.append(line[0])
            elm_type.append(line[1]) 
            if gmshV == 4:
                phys_entity.append(phys)
            else:
                phys_entity.append(line[4]) 
            
            for j in range(npere):
                con_matrix[j].append(line[st+j]-1)
        else:
            ignored_elements += 1
            
    print("ignoring %i elements in the mesh file, as they are not required for R2/R3t"%ignored_elements)
    
    real_no_elements=len(nat_elm_num) #'real' number of elements that we actaully want
    if real_no_elements == 0:
        print("no elements found... aborting" )
        raise Exception ("No elements have been imported, please check formatting of .msh file")
        
    sorted_con_matrix = []
    for c in range(npere):
        sorted_con_matrix.append([])
        sorted_con_matrix[c] = [con_matrix[c][i] for i in range(real_no_elements)]
    
    num_corrected = 0
    #make sure nodes are counterclockwise as this is what R2/R3t expects
    if npere ==3: 
        for i in range(real_no_elements):
            x = [nodex[con_matrix[k][i]] for k in range(3)]
            y = [nodey[con_matrix[k][i]] for k in range(3)]
            z = [nodez[con_matrix[k][i]] for k in range(3)]
            if ccw([x[0],y[0]],[x[1],y[1]],[x[2],y[2]]) == 1: # connection matrix needs fiddling 
                sorted_con_matrix[1][i] =  con_matrix[2][i]
                sorted_con_matrix[2][i] =  con_matrix[1][i]
                num_corrected+=1

    elif npere == 4:
        for i in range(real_no_elements):
            x = [nodex[con_matrix[k][i]] for k in range(4)]
            y = [nodey[con_matrix[k][i]] for k in range(4)]
            z = [nodez[con_matrix[k][i]] for k in range(4)]
            if check_tetra(x, y, z) == 1: # connection matrix needs fiddling 
                sorted_con_matrix[1][i] =  con_matrix[2][i]
                sorted_con_matrix[2][i] =  con_matrix[1][i]
                num_corrected+=1

    elif npere == 6: 
        for i in range(real_no_elements):
            x = [nodex[con_matrix[k][i]] for k in range(6)]
            y = [nodey[con_matrix[k][i]] for k in range(6)]
            z = [nodez[con_matrix[k][i]] for k in range(6)]
            #see if triangle is counter-clockwise
            
            if ccw([x[0],y[0]],[x[1],y[1]],[x[2],y[2]]) == 1: #points are clockwise and therefore need swapping round
                sorted_con_matrix[1][i] =  con_matrix[2][i]
                sorted_con_matrix[2][i] =  con_matrix[1][i]
                sorted_con_matrix[4][i] =  con_matrix[5][i]
                sorted_con_matrix[5][i] =  con_matrix[4][i]
                num_corrected=+1
            
    
    print('%i Node orderings had to be corrected into a counterclockwise direction'%num_corrected) 
            
    elm_id = [i+1 for i in range(real_no_elements)]
            
    mesh_dict = {'num_elms':real_no_elements,
                'num_nodes':no_nodes,
                'dump':dump,      
                'node_x':nodex,#x coordinates of nodes 
                'node_y':nodey,#y coordinates of nodes
                'node_z':nodez,#z coordinates of nodes 
                'node_id':node_num,#node id number 
                'elm_id':elm_id,#element id number 
                'node_data':sorted_con_matrix,#nodes of element vertices
                'parameters':phys_entity,#the values of the attributes given to each cell 
                'parameter_title':'material',
                'dict_type':'mesh_info',
                'original_file_path':file_path} 
    
    print('Finished reading .msh file')
    
    return mesh_dict # return a python dictionary 

#write a dat file 
def write_dat(mesh_dict,file_path='mesh.dat', param=None, zone=None):
    """ Write a mesh.dat kind of file for mesh input for R2. R2 takes a mesh
    input file for triangle meshes, so this function is only relevant for
    triangle meshes.
    
    Parameters
    ------------
    file_path: string, optional
        Path to the file. By default 'mesh.dat' is saved in the working directory. 
    zone: array like, optional
        An array of integers which are assocaited with regions/materials in the mesh. 
        Useful in the case of an inversion which has a boundary constraint. 
    param: array-like, optional
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
    
    con_matrix = mesh_dict['node_data']
    npere = len(con_matrix)
    num_elms = mesh_dict['num_elms']
    num_nodes = mesh_dict['num_nodes']
    con_matrix = mesh_dict['node_data']
    elm_id = mesh_dict['elm_id']
    
    #write to mesh.dat total num of elements and nodes
    if npere!=3:
        fid.write('%i %i 1 0 %i\n'%(num_elms,num_nodes,npere))
    else:
        fid.write('%i %i 0\n'%(num_elms,num_nodes))

    #compute zones if present 
    if zone  is None:
        zone = [1]*num_elms# default zone = 1 
    else:
        if len(zone) != num_elms:
            raise IndexError("the number of zone parameters does not match the number of elements")
        elif min(zone) == 0:
            zone = [zone[i] +1 for i in range(num_elms)]
    
    if param  is None:
        param = [i+1 for i in range(num_elms)]
    else:
        if len(param) != num_elms:
            raise IndexError("the number of parameters does not match the number of elements")
    
    #write out elements         
    for i in range(num_elms):
        fid.write("%i "%elm_id[i])
        [fid.write("%i "%(con_matrix[k][i]+1)) for k in range(npere)]
        fid.write("%i %i\n"%(param[i],zone[i]))

    #now add nodes
    nodex = mesh_dict['node_x']
    nodey = mesh_dict['node_y']
    nodez = mesh_dict['node_z']
    if npere!=3:
        for i in range(num_nodes):
            ni_no=i+1
            fid.write("%i %6.3f %6.3f %6.3f\n"%#node number, x coordinate, y coordinate, z coordinate
                      (ni_no,
                       nodex[i],
                       nodey[i],
                       nodez[i]))
        fid.write('1')
    else:
        for i in range(num_nodes):
            ni_no=i+1
            fid.write("%i %6.3f %6.3f\n"%#node number, x coordinate, z coordinate
                      (ni_no,
                       nodex[i],
                       nodez[i]))

    fid.close()#close the file 
    print('written mesh.dat file to \n%s'%file_path)
    
#%% run script (only if ran as primary module)
if __name__ == "__main__":
    #setup optional command line arguments 
    parser = argparse.ArgumentParser(description=Desc,epilog=note)
    parser.add_argument('--msh',type=str,help='Path to mesh file output from gmsh')
    parser.add_argument('--dat',type=str,help='Path to output .dat file')
    parser.add_argument('--zone',type=bool,help='Make "True" to keep physical entities as mesh zones')
    arg = parser.parse_args()
    
    print("############## .msh to .dat format #############")
    
    if arg.zone is None:
        arg.zone = False # by default ignore the physical entity 
        print('Gmsh Physical Entity will be included as mesh zones')
    else:
        print('Gmsh Physical Entity will be ignored')
    
    
    if arg.msh is None: # no command line argument given  so use UI 
        print("Please select a .msh file to convert...")
        fpath = ui_open_file() # open mesh file 
    else:
        print('(Path to gmsh file given in command line)')
        fpath = arg.msh
    
    mesh = msh_parse(fpath) # parse it 
    print("########### File imported successfully! #########")
    
    if arg.dat is None:
        print('Please choose where to save the exported .dat file')
        fsave = ui_new_file() # choose where to save the output 
    else:
        print('(.dat file path given in command line)')
        fsave = arg.dat   
            
    if fsave.find('mesh3d.dat')==-1 and fsave.find('mesh.dat')==-1:
        print('Note that R2/R3t will not read in the mesh file unless its called "mesh.dat"/"mesh3d.dat"')
        
    if arg.zone:
        write_dat(mesh,file_path=fsave,zone=mesh['parameters'])
    else:
        write_dat(mesh,file_path=fsave) # write .dat file 
    print("############### All okay- END. ##################")
    time.sleep(4)
    
