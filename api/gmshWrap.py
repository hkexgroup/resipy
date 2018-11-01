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
    moving_average() - as says on tin, used to smooth surface topography which is repeated at the base of the fine mesh region in the inversion
    genGeoFile () - generates a .geo file for gmsh
    gmsh2R2mesh () - converts a gmsh.msh file to a mesh.dat file readable by R2
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
#pyR2 library 
#try:
from api.isinpolygon import isinpolygon 
#except ModuleNotFoundError:
#    from isinpolygon import isinpolygon 

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
def genGeoFile(electrodes, electrode_type = None, geom_input = None,
               file_path='mesh.geo',doi=-1,cl=-1,cl_factor=2):
    """
    writes a gmsh .geo file for a 2d study area with topography assuming we wish to add electrode positions
    
    Parameters
    ----------
    geom_input: dict
        a dictionary of electrode coordinates, surface topography, 
                    #borehole electrode coordinates, and boundaries 
    file_path: string, optional 
        name of the generated gmsh file (can include file path also) (optional)
    doi: float, optional 
        depth of investigation (optional) (in meters)
    cl: float, optional
        characteristic length (optional), essentially describes how big the nodes 
        assocaited elements will be. Usually no bigger than 5. 
    cl_factor: float, optional 
        This allows for tuning of the incrimental size increase with depth in the 
        mesh, usually set to 2 such that the elements at the DOI are twice as big as those
        at the surface. The reasoning for this is becuase the sensitivity of ERT drops
        off with depth. 
    
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
            
            geom_input = {'surface': [surf_x,surf_z],
              'electrode':[elec_x,elec_z],
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
#    if 'electrode' not in geom_input:
#        elec_x=[]
#        elec_z=[]
#        topo_x = geom_input['surface'][0]
#        topo_z = geom_input['surface'][1]
    if 'surface' not in geom_input:
        print("surface not in geom_input")
        elec_x = geom_input['electrode'][0]
        elec_z = geom_input['electrode'][1]
        topo_x = [elec_x[0] - 5*np.mean(np.diff(elec_x)),
                  elec_x[-1] + 5*np.mean(np.diff(elec_x))]
        topo_z = [elec_z[0],elec_z[-1]]
    else:
        topo_x = geom_input['surface'][0]
        topo_z = geom_input['surface'][1]       
        elec_x = geom_input['electrode'][0]
        elec_z = geom_input['electrode'][1]
    
        
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
        if bh_flag is True:
            doi = abs(min(geom_input['borehole1'][1])) + abs(0.05*min(geom_input['borehole1'][1]))
        elif bu_flag is True:
            doi = abs(min(geom_input['buried1'][1]))*2# + abs(0.05*min(geom_input['buried1'][1]))
        else:
            doi = abs(np.max(elec_x) - np.min(elec_x))/2
    if cl == -1:
        if bh_flag:
            cl = abs(np.mean(np.diff(geom_input['borehole1'][1]))/2)
        else:
            cl = np.mean(np.diff(elec_x))/2
            
    if len(topo_x) != len(topo_z):
        raise ValueError("topography x and y arrays are not the same length!")
    if len(elec_x) != len(elec_z):
        raise ValueError("electrode x and y arrays are not the same length!")
    
    if file_path.find('.geo')==-1:
        file_path=file_path+'.geo'#add file extension if not specified already
    
    #start to write the file  
    fh = open(file_path,'w') #file handle
    
    fh.write("//Jamyd91's gmsh wrapper code version 0.8 (run the following in gmsh to generate a triangular mesh with topograpghy)\n")
    fh.write("//2D mesh coordinates\n")
    fh.write("cl=%.2f;//define characteristic length\n" %cl)
    fh.write("//Define surface points\n")
    #we have surface topograpghy, and electrode positions to make use of here:
    x_pts=np.append(topo_x,elec_x)#append our electrode positions to our topograpghy 
    y_pts=np.append(topo_z,elec_z)
    flag=['topography point']*len(topo_x)
    flag=flag+(['electrode location']*len(elec_x))   
    #deal with end case electrodes 
    if min(elec_x) == min(x_pts):
        x_pts = np.append(x_pts,elec_x[0] - 5*np.mean(np.diff(elec_x)))
        y_pts = np.append(y_pts,elec_z[0])
        flag.append('topography point')
    if max(elec_x) == max(x_pts):
        x_pts = np.append(x_pts,elec_x[-1] + 5*np.mean(np.diff(elec_x)))
        y_pts = np.append(y_pts,elec_z[-1])
        flag.append('topography point')
     
    idx=np.argsort(x_pts)#now resort the x positions in ascending order
    x_pts=x_pts[idx]#compute sorted arrays
    y_pts=y_pts[idx]
    z_pts=np.zeros((len(x_pts),1))#for now let z = 0 
    #add flags which distinguish what each point is 
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
    sur_pnt_cache=[] # surface points cache 
    for i in range(len(x_pts)):
        tot_pnts=tot_pnts+1
        fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl};//%s\n"%(tot_pnts,x_pts[i],y_pts[i],z_pts[i],flag_sort[i]))
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
    
    #reflect surface topography at base of fine mesh area. 
    x_base = x_pts[::2]
    z_base = moving_average(y_pts[::2] - abs(doi),N=5) # compute the depth to the points at the base of the survey, + downsample
    if len(x_pts)%2 == 0:#bug fix
        z_base = np.append(z_base,y_pts[-1]- abs(doi))#puts in extra point at base underneath last x and y point
        x_base = np.append(x_base,x_pts[-1])
    # a smoothed version of the topography ... 
    
    basal_pnt_cache = []
    for i in range(len(x_base)):
        tot_pnts=tot_pnts+1
        fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl*%.2f};//%s\n"%(tot_pnts,x_base[i],z_base[i],0,cl_factor,
                 'base of smoothed mesh region'))
        basal_pnt_cache.append(tot_pnts)
    
    fh.write("//make a polygon by defining lines between points just made.\n")
    basal_ln_cache=[]
    for i in range(len(x_base)-1):
        tot_lins=tot_lins+1
        fh.write("Line(%i) = {%i,%i};\n"%(tot_lins,basal_pnt_cache[i],basal_pnt_cache[i+1]))
        basal_ln_cache.append(tot_lins)
    
    fh.write("//Add lines at the end points of each of the fine mesh region.\n")
    tot_lins=tot_lins+1;# add to the number of lines rolling total 
    fh.write("Line(%i) = {%i,%i};\n"%(tot_lins,sur_pnt_cache[0],basal_pnt_cache[0]))#line from first point on surface to depth
    tot_lins=tot_lins+1
    fh.write("Line(%i) = {%i,%i};\n"%(tot_lins,sur_pnt_cache[-1],basal_pnt_cache[-1]))#line going bottom to last electrode point
    end_ln_cache=[tot_lins-1,tot_lins]
    
    #compile line numbers into a line loop.
    fh.write("//compile lines into a line loop for a mesh surface/region.\n")
    sur_ln_cache_flipped = list(np.flipud(np.array(sur_ln_cache))*-1)
    fine_msh_loop = [end_ln_cache[0]] + basal_ln_cache + [-1*end_ln_cache[1]] + sur_ln_cache_flipped
    fh.write("Line Loop(1) = {%s};\n"%str(fine_msh_loop).strip('[').strip(']')) # line loop for fine mesh region 
    fh.write("Plane Surface(1) = {1};//Fine mesh region surface\n")
    
    #now extend boundaries beyond flanks of survey area (so generate your Neummon boundary)
    fh.write("\n//Add background region (Neumann boundary) points\n")
    cl_factor2=50#characteristic length multipleier for Nuemon boundary 
    cl2=cl*cl_factor2#assign new cl, this is so mesh elements get larger from the main model
    fh.write("cl2=%.2f;//characteristic length for background region\n" %cl2)
    #Background region propeties, follow rule of thumb that background should extend 100*electrode spacing
    e_spacing=np.mean(np.diff(elec_x))
    flank=e_spacing*100
    b_max_depth=-abs(doi)-flank#background max depth
    #add nuemon boundaries on left hand side
    n_pnt_cache=[0,0,0,0]#cache for the indexes of the neumon boundary points 
    tot_pnts=tot_pnts+1;n_pnt_cache[0]=tot_pnts
    fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl2};//far left upper point\n"%(tot_pnts,x_pts[0]-flank,y_pts[0],z_pts[0]))
    tot_pnts=tot_pnts+1;n_pnt_cache[1]=tot_pnts
    fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl2};//far left lower point\n"%(tot_pnts,x_pts[0]-flank,b_max_depth,z_pts[0]))
    #add nuemon boundary points on right hand side
    tot_pnts=tot_pnts+1;n_pnt_cache[2]=tot_pnts
    fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl2};//far right upper point\n"%(tot_pnts,x_pts[-1]+flank,y_pts[-1],z_pts[-1]))
    tot_pnts=tot_pnts+1;n_pnt_cache[3]=tot_pnts
    fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl2};//far right lower point\n"%(tot_pnts,x_pts[-1]+flank,b_max_depth,z_pts[-1]))
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
    
    fh.write("//Add line loops and plane surfaces to for nuemon region\n")
    #now add background region line loop (cos this be made more efficent?)
    basal_ln_cache_flipped = list(np.flipud(np.array(basal_ln_cache))*-1)
    coarse_msh_loop = n_ln_cache + [end_ln_cache[1]] + basal_ln_cache_flipped + [-1*end_ln_cache[0]]
    fh.write("Line Loop(2) = {%s};\n"%str(coarse_msh_loop).strip('[').strip(']')) # line loop for fine mesh region 
    fh.write("Plane Surface(2) = {2};//Fine mesh region surface\n")
    
    fh.write("\n//Make a physical surface\n")
    fh.write("Physical Surface(1) = {1, 2};\n")
    
    #now we want to return the point values of the electrodes, as gmsh will assign node numbers to points
    #already specified in the .geo file. This will needed for specifying electrode locations in R2.in   
    node_pos=[i+1 for i, j in enumerate(flag_sort) if j == 'electrode location']
    node_pos = np.array(node_pos)
    
    #add borehole vertices and line segments to the survey mesh
    fh.write("\n//Adding boreholes? \n")
    no_lin=tot_lins#+1
    no_pts=tot_pnts#+1
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
            fh.write("//end of borehole strings.\n")
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
            fh.write("//end of buried electrodes.\n")
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
            fh.write("//end of polygons.\n")
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
            fh.write("//end of boundaries.\n")
            print('%i boundary(ies) added to input file'%(count-1))
            break              
                    
    fh.write("\n//j'ai fini!\n")
    fh.close()
    print("writing .geo to file completed, save location:\n%s\n"%os.getcwd())

    return node_pos,file_path

#%% convert gmsh 2d mesh to R2
def gmsh2R2mesh(file_path='ask_to_open',save_path='default',return_mesh=False,poly_data=None):
    """
    Converts a gmsh mesh file into a mesh.dat file needed for R2. 
    
    Parameters
    ----------
    file_path: string
        file path to mesh file. note that a error will occur if the file format is not as expected
    save_path: string
        leave this as default to save the file in the working directory, make this 'ask_to_open' to open a dialogue box, else enter a custom file path.
    return_mesh: boolian, optional
        make this 'yes' to return a mesh dictionary (not class) 
    poly_data: dict, optional
        Each key in the dictionary should be array (2 lists) with 2 columns with the x and y coordinates of the polygon.   
    
    Returns
    ----------
    mesh_info: dict
        dictionary with some 'useful' info about the mesh
    """
    warnings.warn("gmsh2R2mesh will be depreciated in future versions of pyR2.")
    
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
    
    #compute zones if present 
    zone = np.ones(real_no_elements, dtype=int) # default zone = 1 
    if isinstance(poly_data,dict):
        for i, key in enumerate(poly_data):
            poly_x=poly_data[key][0]#polygon x coordinates
            poly_y=poly_data[key][1]#polygon y coordinates
            inside = isinpolygon(np.array(centriod_x),
                                 np.array(centriod_y),
                                 (poly_x,poly_y))
            zone[inside]=i+2 # add 2 because defualt = 1 and indexing starts at 0 (so 0 + 2 = 2)
            
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
                'parameters':zone,#the values of the attributes given to each cell 
                'parameter_title':'material',
                'dict_type':'mesh_info',
                'original_file_path':file_path} 
            #the information here is returned as a mesh dictionary because it much easier to debug 
            
def msh_parse(file_path):
    """
    Converts a gmsh mesh file into a mesh class used in pyR2
    
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
    if file_path=='ask_to_open':#use a dialogue box to open a file
        print("please select the gmsh mesh file you want to convert.\n")
        root=tk.Tk()
        root.withdraw()
        file_path=filedialog.askopenfilename(title='Select mesh file',filetypes=(("mesh files","*.msh"),("all files","*.*")))
    # open file and read in header lines
    print("parsing gmsh mesh...\n")
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
    
    ### return dictionary which can be converted to mesh class ### 
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
            'parameters':phys_entity,#the values of the attributes given to each cell 
            'parameter_title':'material',
            'dict_type':'mesh_info',
            'original_file_path':file_path} 
    
    #we could return a mesh object here, but a dictionary is easier to debug with spyder, 
    #also we'd need to import the mesh class, and its not a good idea to have modules
    #importing each other, as meshTools has a dependency on gmshWrap. 
        
    
#%% test block 
#import parsers as prs     
##import survey 
#elec, df = prs.res2invInputParser(r'C:\Users\jamyd91\Documents\2PhD_projects\Hollin_Hill\Models\res2d_forward_error.dat')
##elec, df = prs.syscalParser(r'C:\Users\jamyd91\Documents\2PhD_projects\R2gui\Data\example_feild_data.txt')
##slope geometry 
#width=170;#width of slope:
#top=100;#top of slope: 100m above datum
#bottom=50;#base of slope 50m ODM;
#ends_addon=40;# how much to add to the width of the slope
#bottom_addon=10;#how much to add to the bottom of the model 
#X=np.array([-ends_addon,0,width,width+ends_addon]);#compile geometry into arrays
#Y=np.array([bottom,bottom,top,top])
#
##%triangular mesh - create geometry and run gmsh wrapper
#geom_input = {'electrode':[elec.T[0],elec.T[1]]}#,
#              #'surface':[X,Y]}
#
#genGeoFile(geom_input)
##%%
#poly2x=np.array([-40,200,200,-40])
#poly2y=np.array([70,70,100,100])
#poly_data = {'region1':[poly2x,poly2y]}
#mesh_dict = gmsh2R2mesh(file_path='ask_to_open',save_path='default',return_mesh=True,poly_data=poly_data)




