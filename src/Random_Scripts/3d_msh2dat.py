""" 
Convert 3D .msh file into a .dat file for R3t and cR3t. Code is adapted from 
the ResIPy package, an all in one python wrapper for R2, cR2, R3t and cR3t. 
The package even includes a graphical user interface. It can be found at 
https://gitlab.com/hkex/pyr2

Author: Jimmy Boyd
Date: 3rd of December 2019 

Note: The code is designed to work with gmsh version 3.06

"""
#import modules (only using standard libs)
import warnings
import tkinter as tk
from tkinter import filedialog

#%% Define relevant functions

def ui_open_file(title = "Please select gmsh mesh file",file_type=["mesh file","msh"]):
    root=tk.Tk()
    root.withdraw()
    file_path=filedialog.askopenfilename(title=title,filetypes=((file_type[0],"*."+file_type[1]),("all files","*.*")))
    return file_path

def ui_new_file(title = "Please save mesh file as",file_type=["mesh file","msh"]):
    root=tk.Tk()
    root.withdraw()
    file_path=filedialog.asksaveasfilename(title=title,filetypes=((file_type[0],"*."+file_type[1]),("all files","*.*")))
    return file_path

def ccw(p,q,r):#code expects points as p=(x,y) and so on ... 
    """
    checks if points in a triangle are ordered counter clockwise. When using R2,
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
    
def type2VertsNo(cell_type):#converts vtk cell types into number of vertices each element has 
    if int(cell_type[0])==5:#then elements are triangles
        return 3
    elif int(cell_type[0])==8 or int(cell_type[0])==9:#elements are quads
        return 4
    elif int(cell_type[0]) == 11: # elements are voxels
        return 8
    elif int(cell_type[0]) == 10:# elements are tetrahedra 
        return 4
    elif int(cell_type[0]) == 13: # elements are 3d wedges 
        return 6
    #add element types as neccessary 
    else:
        print("WARNING: unrecognised cell type")
        return 0

# parse a .msh file
def msh_parse_3d(file_path):
    """
    Converts a 3d gmsh mesh file into a mesh class used in ResIPy
    
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
    print("parsing gmsh mesh...\n")
    fid=open(file_path,'r')# Open text file
    #Idea: Read Mesh format lines $MeshFormat until $Nodes
    dump=fid.readlines()
    fid.close()
    #check the file is a mesh format
    if dump[0].strip() != '$MeshFormat':#removes formating strings, checks if the file is a gmsh file
        raise ImportError("unrecognised file type...")
    mesh_format=dump[1] 
    if mesh_format.strip() != '2.2 0 8':#warn people that the code was developed with this file format in mind
        print('Warning: the mesh file type version is different to the mesh converter development version ... some errors may occur!\n')   
    
    #find where the nodes start 
    for i, line in enumerate(dump):
        if line.find("$Nodes") != -1:#node flag 
            no_nodes=int(dump[i+1])
            node_start = i+2
        if line.find("$EndNodes") != -1:
            node_end = i
            break # stop the loop, should find the nodes start before the end 
    print('importing node coordinates...')
    #read in number of nodes - at line 5
    #allocate lists for node numbers and coordinates
    node_num=[0]*no_nodes
    x_coord=[0]*no_nodes
    y_coord=[0]*no_nodes
    z_coord=[0]*no_nodes
    #read in node information
    node_idx = 1
    for i in range(node_start, node_end):
        line_info=dump[i].split()
        #convert string info into floats
        line_data=[float(k) for k in line_info]
        previous_node = node_idx 
        node_idx = int(line_data[0])
        if not node_idx < previous_node: # else then its some weird flag used by gmsh, ignore
            node_num[node_idx-1]=node_idx
            x_coord[node_idx-1]=line_data[1]
            y_coord[node_idx-1]=line_data[2]
            z_coord[node_idx-1]=line_data[3]
    
    #### read in elements 
    # find where the elements start 
    for i, line in enumerate(dump):
        if line.find("$Elements") != -1:#node flag 
            #no_elements=int(dump[i+1])#number of elements
            element_start = i+2
        if line.find("$EndElements") != -1:
            element_end = i
            break # stop the loop, should find the nodes start before the end 
    print('reading connection matrix')
    
    #engage for loop - this time we want to filter out elements which are not tetrahedron
    #... looking at the gmsh docs its elements of type 2 we are after (R2 only needs this information) 
    nat_elm_num = []#native element number to gmsh
    elm_type = []#element type
    number_of_tags = []
    phys_entity = []#defines the physical entity type the element is assocaited with
    elem_entity = []#which plane surface the element is assocaited with
    node1 = []#first node of triangle 
    node2 = []
    node3 = []
    node4 = []#last node of tetrahedra
    node5 = []
    node6 = []
    npere = 4 
    ignored_elements=0#count the number of ignored elements
    for i in range(element_start,element_end):
        line_data=[int(k) for k in dump[i].split()]
        if line_data[1]==4:# then its the right element type!
        #convert string info into floats and cache data
            prism = False
            nat_elm_num.append(line_data[0])
            elm_type.append(line_data[1]) 
            number_of_tags.append(line_data[2]) 
            phys_entity.append(line_data[3]) 
            elem_entity.append(line_data[4]) 
            node1.append(line_data[5]-1) #we have to take 1 off here cos of how python indexes lists and tuples
            node2.append(line_data[6]-1) 
            node3.append(line_data[7]-1)
            node4.append(line_data[8]-1)
        elif line_data[1]==6: # prism mesh 
            prism = True
            nat_elm_num.append(line_data[0])
            elm_type.append(line_data[1]) 
            number_of_tags.append(line_data[2]) 
            phys_entity.append(line_data[3]) 
            elem_entity.append(line_data[4]) 
            node1.append(line_data[5]-1) 
            node2.append(line_data[6]-1) 
            node3.append(line_data[7]-1)
            node4.append(line_data[8]-1)
            node5.append(line_data[9]-1)
            node6.append(line_data[10]-1)
        else:
            ignored_elements += 1
    print("ignoring %i non-tetrahedra (or prism) elements in the mesh file, as they are not required for R3t"%ignored_elements)
    real_no_elements=len(nat_elm_num) #'real' number of elements that we actaully want
    
    #compute element centres 
    centriod_x=[]
    centriod_y=[]
    centriod_z=[]
    volumes=[]
    
    if not prism:
        node_dump = [node1,node2,node3,node4]
        cell_typ = [10]
        for i in range(real_no_elements):
            n1=(x_coord[node1[i]],y_coord[node1[i]],z_coord[node1[i]])#define node coordinates
            n2=(x_coord[node2[i]],y_coord[node2[i]],z_coord[node2[i]])
            n3=(x_coord[node3[i]],y_coord[node3[i]],z_coord[node3[i]])
            n4=(x_coord[node4[i]],y_coord[node4[i]],z_coord[node4[i]])
            centriod_x.append(sum((n1[0],n2[0],n3[0],n4[0]))/npere)
            centriod_y.append(sum((n1[1],n2[1],n3[1],n4[1]))/npere)
            centriod_z.append(sum((n1[2],n2[2],n3[2],n4[2]))/npere)
    else:
        #make sure in nodes in triangle are counterclockwise as this is what R3t expects
        c_triangles_bot=[]#'corrected' triangles 
        c_triangles_top=[]#'corrected' triangles 
        num_corrected=0#number of elements that needed 'correcting'
        volumes=[]
        #upside_down = 0 
        for i in range(real_no_elements):
            n1=(x_coord[node1[i]],y_coord[node1[i]],z_coord[node1[i]])#define node coordinates
            n2=(x_coord[node2[i]],y_coord[node2[i]],z_coord[node2[i]])
            n3=(x_coord[node3[i]],y_coord[node3[i]],z_coord[node3[i]])
            n4=(x_coord[node4[i]],y_coord[node4[i]],z_coord[node4[i]])
            n5=(x_coord[node5[i]],y_coord[node5[i]],z_coord[node5[i]])
            n6=(x_coord[node6[i]],y_coord[node6[i]],z_coord[node6[i]])
            #see if triangle is counter-clockwise
            
            if n4[2] < n1[2]:#see if prisms are upside down, flip the nodes if so 
                if ccw(n1,n2,n3) == 1: #points are clockwise and therefore need swapping round
                    c_triangles_top.append((node2[i],node1[i],node3[i]))
                    num_corrected=num_corrected+1
                else:
                    c_triangles_top.append((node1[i],node2[i],node3[i]))
                    
                if ccw(n4,n5,n6) == 1: 
                    c_triangles_bot.append((node5[i],node4[i],node6[i]))
                    num_corrected=num_corrected+1
                else:
                    c_triangles_bot.append((node4[i],node5[i],node6[i]))
            
            else:
                if ccw(n1,n2,n3) == 1: #points are clockwise and therefore need swapping round
                    c_triangles_bot.append((node2[i],node1[i],node3[i]))
                    num_corrected=num_corrected+1
                else:
                    c_triangles_bot.append((node1[i],node2[i],node3[i]))
                    
                if ccw(n4,n5,n6) == 1: 
                    c_triangles_top.append((node5[i],node4[i],node6[i]))
                    num_corrected=num_corrected+1
                else:
                    c_triangles_top.append((node4[i],node5[i],node6[i]))
            
            #compute volume (for a prism this is 0.5*base*height*width)
            base=(((n1[0]-n2[0])**2) + ((n1[1]-n2[1])**2))**0.5
            mid_pt=((n1[0]+n2[0])/2,(n1[1]+n2[1])/2)
            width=(((mid_pt[0]-n3[0])**2) + ((mid_pt[1]-n3[1])**2))**0.5
            height = abs(n1[2] - n4[2])
            volumes.append(0.5*base*width*height)
            centriod_x.append(sum((n1[0],n2[0],n3[0],n4[0],n5[0],n6[0]))/npere)
            centriod_y.append(sum((n1[1],n2[1],n3[1],n4[1],n5[1],n6[1]))/npere)
            centriod_z.append(sum((n1[2],n2[2],n3[2],n4[2],n5[2],n6[2]))/npere)
            
        #node_dump = [node1,node2,node3,node4,node5,node6]
        cell_typ = [13]
        node_dump=[[],[],[],[],[],[]]
        for i in range(real_no_elements):
            node_dump[0].append(c_triangles_bot[i][0])#node 1
            node_dump[1].append(c_triangles_bot[i][1])#node 2
            node_dump[2].append(c_triangles_bot[i][2])#node 3
            node_dump[3].append(c_triangles_top[i][0])#node 5
            node_dump[4].append(c_triangles_top[i][1])#node 4
            node_dump[5].append(c_triangles_top[i][2])#node 6
            
        print('%i node orderings had to be corrected into a counterclockwise direction'%num_corrected)
        
        if min(volumes) < 0:
            warnings.warn('Elements with negative or no volume found, R3t unlikely to work!!!')
            
    elm_id = [i+1 for i in range(real_no_elements)]
            
    mesh_dict = {'num_elms':real_no_elements,
            'num_nodes':no_nodes,
            'dump':dump,      
            'node_x':x_coord,#x coordinates of nodes 
            'node_y':y_coord,#y coordinates of nodes
            'node_z':z_coord,#z coordinates of nodes 
            'node_id':node_num,#node id number 
            'elm_id':elm_id,#element id number 
            'num_elm_nodes':elm_type[0],#number of points which make an element
            'node_data':node_dump,#nodes of element vertices
            'elm_centre':[centriod_x,centriod_y,centriod_z],#centre of elements (x,y,z)
            'elm_areas':volumes,
            'cell_type':cell_typ,
            'parameters':phys_entity,#the values of the attributes given to each cell 
            'parameter_title':'material',
            'dict_type':'mesh_info',
            'original_file_path':file_path} 
    return mesh_dict 

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
    
    ndims = 3
    num_elms = mesh_dict['num_elms']
    num_nodes = mesh_dict['num_nodes']
    cell_typ = mesh_dict['cell_type']
    con_matrix = mesh_dict['node_data']
    elm_id = mesh_dict['elm_id']
    
    #write to mesh.dat total num of elements and nodes
    if ndims==3:
        fid.write('%i %i 1 0 %i\n'%(num_elms,num_nodes,type2VertsNo(cell_typ)))
    else:
        fid.write('%i %i\n'%(num_elms,num_nodes))

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
    no_verts = type2VertsNo(cell_typ)
    for i in range(num_elms):
        fid.write("%i "%elm_id[i])
        [fid.write("%i "%(con_matrix[k][i]+1)) for k in range(no_verts)]
        fid.write("%i %i\n"%(param[i],zone[i]))

    #now add nodes
    x_coord = mesh_dict['node_x']
    y_coord = mesh_dict['node_y']
    z_coord = mesh_dict['node_z']
    if ndims==3:
        for i in range(num_nodes):
            ni_no=i+1
            fid.write("%i %6.3f %6.3f %6.3f\n"%#node number, x coordinate, y coordinate, z coordinate
                      (ni_no,
                       x_coord[i],
                       y_coord[i],
                       z_coord[i]))
        fid.write('1')
    else:
        for i in range(num_nodes):
            ni_no=i+1
            fid.write("%i %6.3f %6.3f\n"%#node number, x coordinate, y coordinate
                      (ni_no,
                       x_coord[i],
                       z_coord[i]))

    fid.close()#close the file 
    print('written mesh.dat file to \n%s'%file_path)
    
#%% run script (only if ran as primary module)
if __name__ == "__main__":
    fpath = ui_open_file() # open mesh file 
    mesh = msh_parse_3d(fpath) # parse it 
    fsave = ui_new_file() # choose where to save the output 
    write_dat(mesh,file_path=fsave,zone=mesh['parameters']) # write .dat file 