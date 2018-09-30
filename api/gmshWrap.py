# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 14:32:14 2018 in python 3.6.5
Wrapper for creating 2d triangular meshes with gmsh and converting it 
into a mesh.dat format for R2. The program tri_mesh () expects a directory "exe"
to be within the api (or working) directory with a gmsh.exe inside it. 

@author: jimmy Boyd - jamyd91@bgs.ac.uk
Programs:
    arange () - creates a list of values like np.arange does
    ccw() - checks cartesian points are oreintated clockwise 
    genGeoFile () - generates a .geo file for gmsh
    genGeoFile_adv () - more advanced version of genGeoFile which allows for greater flexiblity 
    gmsh2R2mesh () - converts a gmsh.msh file to a mesh.dat file readable by R2
    tri_mesh () - combines GenGeoFile and gmsh2R2msh functions into one function, returns a mesh object 
    isinpolygon () - checks to see if a point lies inside a polygon

Dependencies: 
    numpy (conda library)
    python3 standard libs
"""
#python standard libraries 
#import tkinter as tk
#from tkinter import filedialog
import os, warnings
#general 3rd party libraries
import numpy as np

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
    
# triangle centriod 
def tri_cent(p,q,r):
    """
    Compute the centre coordinates for a 2d triangle given the x,y coordinates 
    of the vertices.
    program code expects points as p=(x,y) and so on (counter clockwise prefered)
    
    Parameters
    ----------
    
    p - tuple, list,
        The x y coordinates of the point/vertex ie. (x,y) in desired order
    q - " 
    r - " 
    
    Returns
    ----------
    (x,y): tuple    
    """
    Xm=(p[0]+q[0])/2
    Ym=(p[1]+q[1])/2
    k=2/3
    Xc=r[0]+(k*(Xm-r[0]))
    Yc=r[1]+(k*(Ym-r[1]))
    return(Xc,Yc)

#%% write a .geo file for reading into gmsh with topography (and electrode locations)
def genGeoFile(topo_x,topo_y,elec_x,elec_y,file_name="default",doi=-1,cl=-1,path='default',
               bh_survey=False):
    """
    DEPRECIATED
#writes a gmsh .geo file for a 2d study area with topography assuming we wish to add electrode positions
#INPUT:
    #topo_x - x coordinates of the surface
    #topo_y - y coordinates of the surface 
    #elec_x - x coordinates of the electrodes 
    #elec_y - y coordinates of the electrodes
    #file_name - name of the generated gmsh file (optional)
    #doi - depth of investigation (optional)
    #cl - characteristic length (optional)
    #path - directory to save the .geo file 
    #bh_flag = False
#OUTPUT:
    #gmsh . geo file which can be run / loaded into gmsh
###############################################################################
    """
    #at some point were looking to change genGeoFile_adv to the main 2d gmsh
    #input generator, so warn the user. 
    warnings.warn("Deprecaition warning: genGeoFile will be changing in the future to require a dictionary as input")
    
    #formalities and error checks
    if doi == -1:
        doi = abs(np.max(elec_x) - np.min(elec_x))/2
        # print(doi)
        # TODO very rough, better to consider 2/3 of the longest dipole
    if cl == -1:
        cl = np.mean(np.diff(elec_x))/2
    if len(topo_x) != len(topo_y):
        raise ValueError("topograpghy x and y arrays are not the same length!")
    if len(elec_x) != len(elec_y):
        raise ValueError("electrode x and y arrays are not the same length!")
    if file_name.find('.geo')==-1:
        file_name=file_name+'.geo'#add file extension if not specified already
    #start to write the file
    if path=='default':#written to working directory 
        fh=open(file_name,'w')
    else:
        fh=open(os.path.join(path,file_name),'w')#file handle
    fh.write("//Jamyd91's gmsh wrapper code version 0.1 (run the following in gmsh to generate a triangular mesh with topograpghy)\n")
    fh.write("//2D mesh coordinates\n")
    fh.write("cl=%.2f;//define characteristic length\n" %cl)
    fh.write("//Define surface points\n")
    #we have surface topograpghy, and electrode positions to make use of here:
    x_pts=np.append(topo_x,elec_x)#append our electrode positions to our topograpghy 
    y_pts=np.append(topo_y,elec_y)
    idx=np.argsort(x_pts)#now resort the x positions in ascending order
    x_pts=x_pts[idx]#compute sorted arrays
    y_pts=y_pts[idx]
    z_pts=np.zeros((len(x_pts),1))#for now let z = 0 
    #add flags which distinguish what each point is 
    flag=['topography point']*len(topo_x)
    flag=flag+(['electrode location']*len(elec_x))
    flag_sort=[flag[i] for i in idx]
    #we need protection against repeated points, as this will throw up an error in R2 when it comes calculating element areas
    cache_idx=[]
    for i in range(len(x_pts)-1):
        if x_pts[i]==x_pts[i+1] and y_pts[i]==y_pts[i+1]:
            cache_idx.append(i)
    #if duplicated points were dectected we should remove them
    if len(cache_idx)>0:
        x_pts=np.delete(x_pts,cache_idx)#deletes first instance of the duplicate       
        y_pts=np.delete(y_pts,cache_idx)
        print("%i duplicated coordinate(s) were deleted" %len(cache_idx))
        for k in range(len(cache_idx)):#what are we actually deleting?
            if flag_sort[cache_idx[k]] == 'electrode location':
                flag_sort[cache_idx[k]+1] = 'electrode location'
                #this overwrites the flag string to say that this point is an electrode
                #otherwise we'll get a mismatch between the number of electrodes and mesh nodes assigned to the electrodes
        flag_sort=np.delete(flag_sort,cache_idx).tolist()
    #now add the surface points to the file
    tot_pnts=0#setup a rolling total
    for i in range(len(x_pts)):
        val=i+1
        fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl};//%s\n"%(val,x_pts[i],y_pts[i],z_pts[i],flag_sort[i]))
        tot_pnts=tot_pnts+1
    #make the lines between each point
    fh.write("//construct lines between each surface point\n")
    tot_lins=0
    for i in range(len(x_pts)-1):
        val=i+1
        fh.write("Line(%i) = {%i,%i};\n"%(val,val,val+1))
        tot_lins=tot_lins+1
    fh.write("//add points below surface to make a polygon\n")#okay so we want to add in the lines which make up the base of the slope
    max_depth=min(y_pts)-doi
    fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl*1.5};\n"%(tot_pnts+1,x_pts[0],max_depth,z_pts[0]))#point below left hand side of sudy area
    fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl*1.5};\n"%(tot_pnts+2,x_pts[-1],max_depth,z_pts[-1]))#point below right hand side of sudy area
    fh.write("//make a polygon by defining lines between points just made\n")
    fh.write("Line(%i) = {%i,%i};\n"%(tot_lins+1,1,tot_pnts+1))#line from first point on surface to depth
    fh.write("Line(%i) = {%i,%i};\n"%(tot_lins+2,tot_pnts+1,tot_pnts+2))#line going from the 2 new points
    fh.write("Line(%i) = {%i,%i};\n"%(tot_lins+3,tot_pnts+2,tot_pnts))#line going bottom to last electrode point
    #now extend boundaries beyond flanks of slope(so generate your Neummon boundary)
    fh.write("//Add background region (Neumann boundary) points\n")
    cl_factor=150#characteristic length multipleier for Nuemon boundary 
    cl2=cl*cl_factor#assign new cl, this is so mesh elements get larger from the main model
    fh.write("cl2=%.2f;//characteristic length for background region\n" %cl2)
    #Background region propeties, follow rule of thumb that background should extend 100*electrode spacing
    e_spacing=np.mean(np.diff(elec_x))
    flank=e_spacing*100
    b_max_depth=max_depth-flank#background max depth
    #add nuemon boundaries on left hand side
    fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl2};\n"%(tot_pnts+3,x_pts[0]-flank,y_pts[0],z_pts[0]))
    fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl2};\n"%(tot_pnts+4,x_pts[0]-flank,b_max_depth,z_pts[0]))
    #add nuemon boundary points on right hand side
    fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl2};\n"%(tot_pnts+6,x_pts[-1]+flank,y_pts[-1],z_pts[-1]))
    fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl2};\n"%(tot_pnts+5,x_pts[-1]+flank,b_max_depth,z_pts[-1]))
    #make lines encompassing all the points - counter clock wise fashion
    fh.write("//make lines encompassing all the background points - counter clock wise fashion\n")
    fh.write("Line(%i) = {%i,%i};\n"%(tot_lins+4,1,tot_pnts+3))
    fh.write("Line(%i) = {%i,%i};\n"%(tot_lins+5,tot_pnts+3,tot_pnts+4))
    fh.write("Line(%i) = {%i,%i};\n"%(tot_lins+6,tot_pnts+4,tot_pnts+5))
    fh.write("Line(%i) = {%i,%i};\n"%(tot_lins+7,tot_pnts+5,tot_pnts+6))
    fh.write("Line(%i) = {%i,%i};\n"%(tot_lins+8,tot_pnts+6,tot_pnts))
    fh.write("//Add line loops and plane surfaces to mesh\n")
    #make first line loop refering to the model region
    fh.write("Line Loop(%i) = {"%(tot_lins+9))
    for i in np.arange(1,tot_lins,1):
        val=i+1
        fh.write("%i, "%val)
    fh.write("%i, %i, %i, 1};\n"%(-(tot_lins+3),-(tot_lins+2),-(tot_lins+1)))
        #now add background region line loop (cos this be made more efficent?)
    fh.write("Line Loop(%i) = {%i ,%i ,%i ,%i ,%i ,%i ,%i ,%i};\n"%(tot_lins+10,
             tot_lins+1,
             tot_lins+2,
             tot_lins+3,
             -(tot_lins+8),
             -(tot_lins+7),
             -(tot_lins+6),
             -(tot_lins+5),
             -(tot_lins+4)))
    fh.write("//Assign plane surfaces to be meshed\n")
    fh.write("Plane Surface(1) = {%i};\n"%(tot_lins+9))
    fh.write("Plane Surface(2) = {%i};\n"%(tot_lins+10))
    fh.write("//Finally make a physical surface\n")
    fh.write("Physical Surface(1) = {1, 2};\n")
    fh.write("//j'ai fini!")
    fh.close()
    print("writing .geo to file completed\n")
   #now we want to return the point values of the electrodes, as gmsh will assign node numbers to points
   #already specified in the .geo file. This will needed for specifying electrode locations in R2.in   
    node_pos=[i+1 for i, j in enumerate(flag_sort) if j == 'electrode location']
    return node_pos,file_name

#%% write a .geo file for reading into gmsh with topography (and electrode locations)
def genGeoFile_adv(geom_input,file_name="default",doi=-1,cl=-1,path='default'):
    """
    writes a gmsh .geo file for a 2d study area with topography assuming we wish to add electrode positions
    
    Parameters
    ----------
    geom_input: dict
        a dictionary of electrode coordinates, surface topography, 
                    #borehole electrode coordinates, and boundaries 
    file_name: string
        name of the generated gmsh file (can include file path also) (optional)
    doi: float
        depth of investigation (optional) (in meters)
    cl: float
        characteristic length (optional), essentially describes how big the nodes 
        assocaited elements will be. Usually no bigger than 5
    path: string
        directory to save the .geo file
    
    Returns
    ----------
    gmsh .geo file which can be ran in gmsh

    NOTES
    ----------
     geom_input format:
        the code will cycle through numerically ordered keys (strings referencing objects in a dictionary"),
        currently the code expects a 'surface' and 'electrode' key for surface points and electrodes.
        the first borehole string should be given the key 'borehole1' and so on. The code stops
        searching for more keys when it cant find the next numeric key. Same concept goes for adding boundaries
        and polygons to the mesh. See below example:
            
            geom_input = {'surface': [surf_x,surf_y],
              'electrode':[elec_x,elec_y],
              'borehole1':[string1x,string1y],
              'borehole2':[string2x,string2y],
              'boundary1':[bound1x,bound1y],
              'polygon1':[poly1x,poly1y],
              'buried1:[x_coord,y_coord]} 
            
    The code expects that all polygons, boundaries and electrodes fall within x values 
    of the actaul survey area. So make sure your topography / surface electrode points cover 
    the area you are surveying, otherwise some funky errors will occur in the mesh. 

    #### TODO: search through each set of points and check for repeats 
    #### TODO: check with some real data
    #### TODO: change nuemon boundary distance to be based on the survey extent
    #### TODO: add points around the survey to provide some sort of padding? 
    """
    print('Generating gmsh input file...\n')
    #formalities and error checks
    if not isinstance(geom_input,dict):
        raise TypeError ("'geom_input' is not a dictionary type object. Dict type is expected for the first argument of genGeoFile_adv")
    
    #as it stands the program requires some surface points to start with, either in the surface
    #or electrode keys. 
    if 'surface' not in geom_input and 'electrode' not in geom_input:
        raise Exception("niether surface electrode or elevation points have been given to genGeoFile. Aborting... ")
    
    #determine     
    if 'electrode' not in geom_input:
        elec_x=[]
        elec_y=[]
        topo_x = geom_input['surface'][0]
        topo_y = geom_input['surface'][1]
    elif 'surface' not in geom_input:
        elec_x = geom_input['electrode'][0]
        elec_y = geom_input['electrode'][1]
        topo_x = [elec_x[0] - 5*np.mean(np.diff(elec_x)),
                  elec_x[-1] + 5*np.mean(np.diff(elec_x))]
        topo_y = [elec_y[0],elec_y[-1]]
    else:
        topo_x = geom_input['surface'][0]
        topo_y = geom_input['surface'][1]       
        elec_x = geom_input['electrode'][0]
        elec_y = geom_input['electrode'][1]
        
    if 'borehole1' in geom_input: # do we have boreholes ? 
        bh_flag = True
        print('Borehole electrodes present!')
    else: 
        bh_flag = False
        
    if 'buried1' in geom_input: # do we have boreholes ? 
        bu_flag = True
        print('Buried electrodes present!')
    else: 
        bu_flag = False
  
    if doi == -1:#then set to a default 
        if bh_flag:
            doi = abs(min(geom_input['borehole1'][1])) + abs(0.05*min(geom_input['borehole1'][1]))
        else:
            doi = abs(np.max(elec_x) - np.min(elec_x))/2
    if cl == -1:
        if bh_flag:
            cl = abs(np.mean(np.diff(geom_input['borehole1'][1]))/2)
        else:
            cl = np.mean(np.diff(elec_x))/2
            
    if len(topo_x) != len(topo_y):
        raise ValueError("topography x and y arrays are not the same length!")
    if len(elec_x) != len(elec_y):
        raise ValueError("electrode x and y arrays are not the same length!")
    if file_name.find('.geo')==-1:
        file_name=file_name+'.geo'#add file extension if not specified already
    
    #start to write the file
    if path=='default':#written to working directory 
        file_path = file_name
    else:
        file_path = os.path.join(path,file_name)
    
    fh = open(file_path,'w') #file handle
    
    fh.write("//Jamyd91's gmsh wrapper code version 0.2 (run the following in gmsh to generate a triangular mesh with topograpghy)\n")
    fh.write("//2D mesh coordinates\n")
    fh.write("cl=%.2f;//define characteristic length\n" %cl)
    fh.write("//Define surface points\n")
    #we have surface topograpghy, and electrode positions to make use of here:
    x_pts=np.append(topo_x,elec_x)#append our electrode positions to our topograpghy 
    y_pts=np.append(topo_y,elec_y)
    idx=np.argsort(x_pts)#now resort the x positions in ascending order
    x_pts=x_pts[idx]#compute sorted arrays
    y_pts=y_pts[idx]
    z_pts=np.zeros((len(x_pts),1))#for now let z = 0 
    #add flags which distinguish what each point is 
    flag=['topography point']*len(topo_x)
    flag=flag+(['electrode location']*len(elec_x))
    flag_sort=[flag[i] for i in idx]
    
    #we need protection against repeated points, as this will throw up an error in R2 when it comes calculating element areas
    cache_idx=[]
    for i in range(len(x_pts)-1):
        if x_pts[i]==x_pts[i+1] and y_pts[i]==y_pts[i+1]:
            cache_idx.append(i)
    
    #if duplicated points were dectected we should remove them
    if len(cache_idx)>0:
        x_pts=np.delete(x_pts,cache_idx)#deletes first instance of the duplicate       
        y_pts=np.delete(y_pts,cache_idx)
        print("%i duplicated coordinate(s) were deleted" %len(cache_idx))
        for k in range(len(cache_idx)):#what are we actually deleting?
            if flag_sort[cache_idx[k]] == 'electrode location':
                flag_sort[cache_idx[k]+1] = 'electrode location'
                #this overwrites the flag string to say that this point is an electrode
                #otherwise we'll get a mismatch between the number of electrodes and mesh nodes assigned to the electrodes
        flag_sort=np.delete(flag_sort,cache_idx).tolist()
    
    #now add the surface points to the file
    print('adding surface points and electrodes to input file...')
    tot_pnts=0#setup a rolling total for points numbering
    for i in range(len(x_pts)):
        val=i+1
        fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl};//%s\n"%(val,x_pts[i],y_pts[i],z_pts[i],flag_sort[i]))
        tot_pnts=tot_pnts+1
    
    #make the lines between each point
    fh.write("//construct lines between each surface point\n")
    tot_lins=0
    for i in range(len(x_pts)-1):
        val=i+1
        fh.write("Line(%i) = {%i,%i};\n"%(val,val,val+1))
        tot_lins=tot_lins+1
    
    fh.write("//add points below surface to make a polygon\n")#okay so we want to add in the lines which make up the base of the slope
    max_depth=-doi
    cl_factor = 2  
    if bh_flag:
        cl_factor = 1       
    fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl*%f};\n"%(tot_pnts+1,x_pts[0],max_depth,z_pts[0],cl_factor))#point below left hand side of sudy area
    fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl*%f};\n"%(tot_pnts+2,x_pts[-1],max_depth,z_pts[-1],cl_factor))#point below right hand side of sudy area
    fh.write("//make a polygon by defining lines between points just made\n")
    fh.write("Line(%i) = {%i,%i};\n"%(tot_lins+1,1,tot_pnts+1))#line from first point on surface to depth
    fh.write("Line(%i) = {%i,%i};\n"%(tot_lins+2,tot_pnts+1,tot_pnts+2))#line going from the 2 new points
    fh.write("Line(%i) = {%i,%i};\n"%(tot_lins+3,tot_pnts+2,tot_pnts))#line going bottom to last electrode point
    #now extend boundaries beyond flanks of survey area (so generate your Neummon boundary)
    fh.write("//Add background region (Neumann boundary) points\n")
    cl_factor=50#characteristic length multipleier for Nuemon boundary 
    cl2=cl*cl_factor#assign new cl, this is so mesh elements get larger from the main model
    fh.write("cl2=%.2f;//characteristic length for background region\n" %cl2)
    #Background region propeties, follow rule of thumb that background should extend 100*electrode spacing
    e_spacing=np.mean(np.diff(elec_x))
    flank=e_spacing*100
    b_max_depth=max_depth-flank#background max depth
    #add nuemon boundaries on left hand side
    fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl2};\n"%(tot_pnts+3,x_pts[0]-flank,y_pts[0],z_pts[0]))
    fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl2};\n"%(tot_pnts+4,x_pts[0]-flank,b_max_depth,z_pts[0]))
    #add nuemon boundary points on right hand side
    fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl2};\n"%(tot_pnts+6,x_pts[-1]+flank,y_pts[-1],z_pts[-1]))
    fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl2};\n"%(tot_pnts+5,x_pts[-1]+flank,b_max_depth,z_pts[-1]))
    #make lines encompassing all the points - counter clock wise fashion
    fh.write("//make lines encompassing all the background points - counter clock wise fashion\n")
    fh.write("Line(%i) = {%i,%i};\n"%(tot_lins+4,1,tot_pnts+3))
    fh.write("Line(%i) = {%i,%i};\n"%(tot_lins+5,tot_pnts+3,tot_pnts+4))
    fh.write("Line(%i) = {%i,%i};\n"%(tot_lins+6,tot_pnts+4,tot_pnts+5))
    fh.write("Line(%i) = {%i,%i};\n"%(tot_lins+7,tot_pnts+5,tot_pnts+6))
    fh.write("Line(%i) = {%i,%i};\n"%(tot_lins+8,tot_pnts+6,tot_pnts))
    fh.write("//Add line loops and plane surfaces to mesh\n")
    #make first line loop refering to the model region
    fh.write("Line Loop(%i) = {"%(tot_lins+9))
    for i in np.arange(1,tot_lins,1):
        val=i+1
        fh.write("%i, "%val)
    fh.write("%i, %i, %i, 1};\n"%(-(tot_lins+3),-(tot_lins+2),-(tot_lins+1)))
        #now add background region line loop (cos this be made more efficent?)
    fh.write("Line Loop(%i) = {%i ,%i ,%i ,%i ,%i ,%i ,%i ,%i};\n"%(tot_lins+10,
             tot_lins+1,
             tot_lins+2,
             tot_lins+3,
             -(tot_lins+8),
             -(tot_lins+7),
             -(tot_lins+6),
             -(tot_lins+5),
             -(tot_lins+4)))
    fh.write("//Assign plane surfaces to be meshed\n")
    fh.write("Plane Surface(1) = {%i};\n"%(tot_lins+9))
    fh.write("Plane Surface(2) = {%i};\n"%(tot_lins+10))
    fh.write("//Make a physical surface\n")
    fh.write("Physical Surface(1) = {1, 2};\n")
    
    #now we want to return the point values of the electrodes, as gmsh will assign node numbers to points
    #already specified in the .geo file. This will needed for specifying electrode locations in R2.in   
    node_pos=[i+1 for i, j in enumerate(flag_sort) if j == 'electrode location']
    node_pos = np.array(node_pos)
    
    #add borehole vertices and line segments to the survey mesh
    fh.write("\n//Adding boreholes? \n")
    no_lin=tot_lins+10
    no_pts=tot_pnts+6
    #add borehole electrode strings
    print('probing for boundaries and other additions to the mesh')
    count = 0
    while True:
        count += 1
        key = 'borehole'+str(count)#count through the borehole keys
        try:
            bhx = geom_input[key][0]#get borehole coordinate information
            bhy = geom_input[key][1]
            try:
                bhz = geom_input[key][2]#allows for future 3d capacity 
            except IndexError:
                bhz = [0]*len(bhx)
            e_pt_idx = [0] *len(bhx)
            fh.write("// string electrodes for borehole %i\n"%count)
            for k in range(len(bhx)):
                no_pts += 1 
                fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl};//borehole %i electrode\n"%(no_pts,bhx[k],bhy[k],bhz[k],count))
                e_pt_idx[k] = no_pts
            fh.write("//put lines between each electrode\n")
            line_idx = []
            node_pos = np.append(node_pos,e_pt_idx) #add borehole nodes to electrode node positions 
            for i in range(len(e_pt_idx)-1):
                idx = e_pt_idx[i]
                no_lin += 1
                fh.write("Line (%i) = {%i,%i};//borehole %i segment\n"%(no_lin,idx,idx+1,count))
                line_idx.append(no_lin)
            fh.write("Line{%s} In Surface{1};\n"%str(line_idx).strip('[').strip(']'))
                
        except KeyError:#run out of borehole keys 
            fh.write("//no more borehole strings to add.\n")
            print('%i boreholes added to input file'%(count-1))
            break
    
    fh.write("\n//Adding buried electrodes? \n")  
    count = 0
    while True:
        count += 1
        key = 'buried'+str(count)#count through the buried keys
        try:
            bhx = geom_input[key][0]#get buried electrode coordinate information
            bhy = geom_input[key][1]
            try:
                bhz = geom_input[key][2]#allows for future 3d capacity 
            except IndexError:
                bhz = [0]*len(bhx)
            e_pt_idx = [0] *len(bhx)
            fh.write("// string electrodes for buried string %i\n"%count)
            for k in range(len(bhx)):
                no_pts += 1 
                fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl};//buried electrode %i\n"%(no_pts,bhx[k],bhy[k],bhz[k],k+1))
                e_pt_idx[k] = no_pts
            node_pos = np.append(node_pos,e_pt_idx) #add borehole nodes to electrode node positions 
            fh.write("Point{%s} In Surface{1};\n"%str(e_pt_idx).strip('[').strip(']'))
                
        except KeyError:#run out of keys 
            fh.write("//no more buried electrode strings to add.\n")
            if bu_flag:
                print('buried electrodes added to input file')
            break
        
    no_plane = 2 # number of plane surfaces so far
    fh.write("\n//Adding polygons?\n")
    count = 0    
    while True:  
        count += 1
        key = 'polygon'+str(count)

        try:
            plyx = geom_input[key][0]
            plyy = geom_input[key][1]
            try:
                plyz = geom_input[key][2]
            except IndexError:
                plyz = [0]*len(plyx)
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
            
        except KeyError:
            fh.write("//no more polygons to add.\n")
            print('%i polygons added to input file'%(count-1))
            break  

    fh.write("\n//Adding boundaries?\n")
    count = 0   
    while True:
        count += 1        
        key = 'boundary'+str(count)

        try:
            bdx = geom_input[key][0]
            bdy = geom_input[key][1]
            try:
                bdz = geom_input[key][2]
            except IndexError:
                bdz = [0]*len(bdx)
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
                
        except KeyError:
            fh.write("//no more boundaries to add.\n")
            print('%i boundary(ies) added to input file'%(count-1))
            break              
                    
    fh.write("\n//j'ai fini!\n")
    fh.close()
    print("writing .geo to file completed, save location:\n%s\n"%file_path)

    return node_pos,file_name

#%% convert gmsh 2d mesh to R2
def gmsh2R2mesh(file_path='ask_to_open',save_path='default',return_mesh=False):
    """
    Converts a gmsh mesh file into a mesh.dat file needed for R2. 
    
    Parameters
    ----------
    file_path: string
        file path to mesh file. note that a error will occur if the file format is not as expected
    save_path: string
        leave this as default to save the file in the working directory, make this 'ask_to_open' to open a dialogue box, else enter a custom file path.
    return_info:
        make this 'yes' if you want to return information from the mesh conversion
    
    Returns
    ----------
    mesh_info: dict
        dictionary with some 'useful' info about the mesh
    """
    if file_path=='ask_to_open':#use a dialogue box to open a file
        print("please select the gmsh mesh file you want to convert.\n")
        root=tk.Tk()
        root.withdraw()
        file_path=filedialog.askopenfilename(title='Select mesh file',filetypes=(("mesh files","*.msh"),("all files","*.*")))
    # open file and read in header lines
    print("converting gmsh mesh into R2 mesh...\n")
    fid=open(file_path,'r')# Open text file
    #Idea: Read Mesh format lines $MeshFormat until $Nodes
    line1=fid.readline()
    #check the file is a mesh format
    if line1.strip() != '$MeshFormat':#removes formating strings, checks if the file is a gmsh file
        raise ImportError("unrecognised file type...")
    mesh_format=fid.readline()#reads the next line
    if mesh_format.strip() != '2.2 0 8':#warn people that the code was developed with this file format in mind
        print('Warning: the mesh file type version is different to the mesh converter development version ... some errors may occur!\n')   
    line3=fid.readline()#endofmeshformat
    line4=fid.readline()#nodes
    
    print('importing node coordinates...')
    #read in number of nodes - at line 5
    no_nodes=int(fid.readline().strip())
    #allocate lists for node numbers and coordinates
    node_num=[0]*no_nodes
    x_coord=[0]*no_nodes
    y_coord=[0]*no_nodes
    z_coord=[0]*no_nodes
    #read in node information
    for i in range(no_nodes):
        line_info=fid.readline().split()
        #convert string info into floats
        data_dump=[float(k) for k in line_info]
        node_num[i]=int(data_dump[0])
        x_coord[i]=data_dump[1]
        y_coord[i]=data_dump[2]
        z_coord[i]=data_dump[3]
    
    #### read in elements   
    print('reading connection matrix')
    #read in two lines $EndNodes and $Elements
    Endnodes=fid.readline()
    Elements=fid.readline()
    #number of elements
    no_elements=int(fid.readline().strip())
    #engage for loop - this time we want to filter out elements which are not triangles
    #... looking at the gmsh docs its elements of type 2 we are after (R2 only needs this information) 
    nat_elm_num = []#native element number to gmsh
    elm_type = []#element type
    number_of_tags = []
    phys_entity = []#defines the physical entity type the element is assocaited with
    elem_entity = []#which plane surface the element is assocaited with
    node1 = []#first node of triangle 
    node2 = []
    node3 = []#last node of triangle 
    ignored_elements=0#count the number of ignored elements
    for i in range(no_elements):
        line_info=fid.readline().split()
        if line_info[1]=='2':# then its the right element type!
        #convert string info into floats and cache data
            data_dump=[int(k) for k in line_info]
            nat_elm_num.append(data_dump[0])
            elm_type.append(data_dump[1]) 
            number_of_tags.append(data_dump[2]) 
            phys_entity.append(data_dump[3]) 
            elem_entity.append(data_dump[4]) 
            node1.append(data_dump[5]) 
            node2.append(data_dump[6]) 
            node3.append(data_dump[7])
        else:
            ignored_elements += 1
    print("ignoring %i non-triangle elements in the mesh file, as they are not required for R2\n"%ignored_elements)
    real_no_elements=len(nat_elm_num) #'real' number of elements that we actaully want
    
    ##clock wise correction and area / centre computations 
    #make sure in nodes in triangle are counterclockwise as this is waht r2 expects
    c_triangles=[]#'corrected' triangles 
    num_corrected=0#number of elements that needed 'correcting'
    centriod_x=[]
    centriod_y=[]
    areas=[]
    for i in range(real_no_elements):
        n1=(x_coord[node1[i]-1],y_coord[node1[i]-1])#define node coordinates
        n2=(x_coord[node2[i]-1],y_coord[node2[i]-1])#we have to take 1 off here cos of how python indexes lists and tuples
        n3=(x_coord[node3[i]-1],y_coord[node3[i]-1])
        #see if triangle is counter-clockwise
        if ccw(n1,n2,n3) == 1: #points are clockwise and therefore need swapping round
            #exchange elements in rows 6 and 7 to change direction
            c_triangles.append((node2[i],node1[i],node3[i]))
            num_corrected=num_corrected+1
        else:
            c_triangles.append((node1[i],node2[i],node3[i]))
        #compute triangle centre
        xy_tuple=tri_cent(n1,n2,n3)#actual calculation
        centriod_x.append(xy_tuple[0])
        centriod_y.append(xy_tuple[1])
        #compute area (for a triangle this is 0.5*base*height)
        base=(((n1[0]-n2[0])**2) + ((n1[1]-n2[1])**2))**0.5
        mid_pt=((n1[0]+n2[0])/2,(n1[1]+n2[1])/2)
        height=(((mid_pt[0]-n3[0])**2) + ((mid_pt[1]-n3[1])**2))**0.5
        areas.append(0.5*base*height)
        
    #print warning if areas of zero found, this will cuase problems in R2
    if min(areas)==0:
        warnings.warn("elements with no area have been detected in 'mesh.dat', inversion with R2 unlikey to work!" )
            
    print("%i element node orderings had to be corrected becuase they were found to be orientated clockwise\n"%num_corrected)
    fid.close()
    ##### write data to mesh.dat kind of file
    if save_path=='default':
        save_path='mesh.dat'
    elif save_path=='ask_to_open':
        print("please select a save location for the converted mesh\n")
        root=tk.Tk()
        root.withdraw()
        save_path=filedialog.asksaveasfilename(title='Select save path',filetypes=(("data files","*.dat"),("all files","*.*")))
    else:
        save_path = os.path.join(save_path,'mesh.dat')
    #open mesh.dat for input      
    fid=open(save_path, 'w')
    print('saving R2 mesh file to:\n%s'%save_path)
    #write to mesh.dat total num of elements and nodes
    fid.write('%i %i\n'%(real_no_elements,no_nodes))
    zone=[1]*real_no_elements
    #add element data following the R2 format
    for i in range(real_no_elements):
        elm_no=i+1
        fid.write("%i %i %i %i %i %i\n"%#element number, nd1, nd2, nd3, parameter,zone.
                  (elm_no,
                   c_triangles[i][0],#node 1
                   c_triangles[i][1],#node 2
                   c_triangles[i][2],#node 3
                   elm_no,#assigning the parameter number as the elm number allows for a unique parameter to be assigned
                   zone[i]))
    #now add nodes
    for i in range(no_nodes):
        ni_no=i+1
        fid.write("%i %6.3f %6.3f\n"%#node number, x coordinate, y coordinate
                  (ni_no,
                   x_coord[i],
                   y_coord[i]))
    fid.close()#close the file 
    print('finished converting mesh!')
    #### find and return some useful information
    if return_mesh:
        print('return mesh info option selected.')
        no_regions=max(elem_entity)#number of regions in the mesh
        regions=arange(1,1,no_regions,1)
        assctns=[]
        #following for loop finds the element number ranges assocaited with a distinct region in the mesh
        for k in regions:
            indx=[m for m in range(len(elem_entity)) if elem_entity[m]==k]
            if len(indx) > 0:
                assctns.append((k,min(indx)+1,max(indx)+1))
        #create a dump of the mesh data incase the user wants to see it later on   
        dump={'nat_elm_num':nat_elm_num,
              'elm_type':elm_type,
              'number_of_tags':number_of_tags,
              'phys_entity':phys_entity,
              'elem_entity':elem_entity,
              'string_data':[line1,mesh_format,line3,line4,Endnodes,Elements]} 
        #convert c_triangles into mesh object format for later recall
        node_dump=[[],[],[]]
        for i in range(real_no_elements):
            node_dump[0].append(c_triangles[i][0]-1)#node 1
            node_dump[1].append(c_triangles[i][1]-1)#node 2
            node_dump[2].append(c_triangles[i][2]-1)#node 3
        #return a dictionary detailing the mesh 
        return {'num_elms':real_no_elements,
                'num_nodes':no_nodes,
                'num_regions':no_regions,
                'element_ranges':assctns,
                'dump':dump,      
                'node_x':x_coord,#x coordinates of nodes 
                'node_y':y_coord,#y coordinates of nodes
                'node_z':z_coord,#z coordinates of nodes 
                'node_id':node_num,#node id number 
                'elm_id':np.arange(1,real_no_elements,1),#element id number 
                'num_elm_nodes':3,#number of points which make an element
                'node_data':node_dump,#nodes of element vertices
                'elm_centre':(centriod_x,centriod_y),#centre of elements (x,y)
                'elm_area':areas,
                'cell_type':[5],
                'parameters':elem_entity,#the values of the attributes given to each cell 
                'parameter_title':'material',
                'dict_type':'mesh_info',
                'original_file_path':file_path} 
            #the information here is returned as a mesh dictionary because it much easier to debug 


#%% gmsh wrapper    - moved to MeshTools    

#%% is in polygon code - originally developed as part of jamyd's slope stability code project. 
def __LineEq(x,y):
    #finds the equation of a line using an inverse operation given an array of
    # x and y data points
    if len(x)==1:#few if statements to check the input is as expected
        raise ValueError ("cant find the equation of a line from a single point!")
    elif len(x)!=len(y):
        raise ValueError ("x and y inputs must be of the same length")
    #now to find the equation of the line
    if sum(np.diff(x))==0:
        #there can be no y intercept the equation of the line takes the form X=constant
        #print("infinite slope detected!")
        return 0,0
    else:
        mat=np.matrix(np.ones((len(x),1)))#vector of ones used in the inverse operator
        x=np.reshape(np.matrix(x),(len(x),1))
        y=np.reshape(np.matrix(y),(len(y),1))#make data into row vectors
        #invert to find model parameters
        G=np.concatenate((x,mat),axis=1)#inverse operator
        mod=((np.transpose(G)*G)**-1)*np.transpose(G)*y#inverse operation to find model
        return float(mod[0]),float(mod[1])#out put is in the form "m,c"

#do the lines then intersect?
def __intersect(A,B,C,D,coord=0):
    #check if line AB intersects line CD. Each coordinate should be passed as a tuple ie. (x,y)
    #pass an extra argment in the function coord==1 if you want to find the intersection of 2 points
    #output is in the form "binary-t/f,(x intercept,y intercept)" 
    if ccw(A,B,C) != ccw(A,B,D) and ccw(C,D,A) != ccw(C,D,B):
        dump=1
    else:
        dump=0
# returns binary (dump) if lines overlap
# see link above for full explanation of the above, note that the code cannot deal with the colinear case
# we may also want to find the intersection point of the 2 lines if they do intersect
    if coord==1 and dump==0:
        return dump,('n/a','n/a')
    elif coord==0:
        return dump
    elif coord==1 and dump==1:
        x_a=(A[0],B[0])#define x and y coordinates 
        y_a=(A[1],B[1])
        x_b=(C[0],D[0])
        y_b=(C[1],D[1])
        #if there is a horizontal and vertical line then we have a special case where we dont need to solve
        #the equations of the lines, we already know the relevant coordinates, also the 
        #LineEq code (above) doesnt deal well with vertical lines
        if np.diff(y_a)==0  and np.diff(x_b)==0:
            return dump,(x_b[0],y_a[0])
        elif np.diff(y_b)==0  and np.diff(x_a)==0:
            return dump,(x_a[0],y_b[0])
        else:
            #if the lines are not both vertical and horizontal we need to solve the 
            #equations of the lines then substitute values in order to find x and y intercept
            m_a,c_a=__LineEq(x_a,y_a)
            m_b,c_b=__LineEq(x_b,y_b)
            x_int=(c_b-c_a)/(m_a-m_b)
            y_int=(m_a*x_int)+c_a
        return dump,(x_int,y_int)  


def isinpolygon(x,y,polydata,*args):
    """
     "Is point inside region?" code
    following solution posted at:
    https://www.geeksforgeeks.org/how-to-check-if-a-given-point-lies-inside-a-polygon/
    
    Parameters
    ----------
    x: float
        x coordinate of query point
    y: float
        y coordinate of query point 
    poly data: list
        a list of polygon data in the form
                #| x coordinates | y coordinates
    *args: optional
        first optional argument changes the distance to which the point is ray casted to the right
    
    Returns
    ----------
    ans: binary 
        1 for true, 0 for false  
    """
    polyx=polydata[0]#polygon x coordinates
    polyy=polydata[1]#polygon y coordinates
    if len(polyx)!=len(polyy):
        raise ValueError('polygon vertices, xy coordinate, arrays are not the same length, check polygon generation scheme!')
    if len(args)==1:
        B=(x+args[0],y)
    else:
        B=(x+10000,y)
    #line AB is the coordinate and raycasted coordinate (B)
    A=(x,y)
    #line CD will represent a segment from the polygon
    count=0#number of times an intersection has been detected 
    for i in range(len(polyx)):
        if i==len(polyx)-1:#its the last iteration
            C=(polyx[-1],polyy[-1])#take last and first polygon point
            D=(polyx[0],polyy[0])
        else:
            C=(polyx[i],polyy[i])
            D=(polyx[i+1],polyy[i+1])
        if __intersect(A,B,C,D)==1:
            count=count+1
    if count==1:
        return True#returns 1 if point is inside polygon 
    elif count==0 or count ==2:
        return False#point must be outside of polygon
    else:
        #this last option shouldnt be possible unless something has gone wrong
        print("something may have gone wrong in the 'isinpolygon' algorithm")
        return 0


#%% test code
#mesh, element_ranges = tri_mesh(np.arange(10), np.zeros(10),
#                np.arange(10), np.zeros(10), keep_files=True, save_path='../test/mesh.dat')
#mesh.show()



