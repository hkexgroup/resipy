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
    msh_parse() - converts a 2d gmsh.msh file to a mesh class readable by R2. 

Dependencies: 
    numpy (conda library)
    python3 standard libs
"""
#python standard libraries 
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
# 2D half space problem 
def genGeoFile(electrodes, electrode_type = None, geom_input = None,
               file_path='mesh.geo',doi=-1,dp_len=-1,cl=-1,cl_factor=2):
    """
    writes a gmsh .geo file for a 2d study area with topography assuming we wish to add electrode positions
    
    Parameters
    ----------
    electrodes: array like
        first column/list is the x coordinates of electrode positions, second column
        is the elevation
    electrode_type: list, optional
        List should be the same length as the electrode coordinate argument. Each entry in
        the list references the type of electrode: 
            'electrode' = surface electrode coordinate, will be used to construct the topography in the mesh
            'buried' = buried electrode, placed the mesh surface
            'borehole#' = borehole electrode, electrodes will be placed in the mesh with a line connecting them. 
                        borehole numbering starts at 1 and ascends numerically by 1. 
            'remote' = remote electrode, this will not be actaully meshed, as often remote electrode coordinates
                        are given arbitary values like -999. Instead every remote electrode found will be 
                        randomly added as a node in the background region towards the lower left hand corner of 
                        of the mesh. 
    geom_input: dict, optional
        Allows for further customisation of the 2D mesh, its a
        dictionary contianing surface topography, polygons and boundaries 
    file_path: string, optional 
        name of the generated gmsh file (can include file path also) (optional)
    doi: float, optional 
        depth of investigation (optional) (in meters)
    cl: float, optional
        characteristic length (optional), essentially describes how big the nodes 
        assocaited elements will be. Usually no bigger than 5. 
    cl_factor: float, optional 
        Characteristic length factor, this allows for tuning of the incrimental 
        size increase with depth in the mesh, usually set to 2 such that the 
        elements at the DOI are twice as big as those at the surface. The reasoning
        for this is because the sensitivity of ERT drops off with depth. 
    dp_len: float, optional 
        Largest dipole length in the 2D array. Default is the largest electrode 
        spacing. Controls the multipier applied to the edges of the coarse mesh region. 
        
    
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

    #### TODO: search through each set of points and check for repeats 
    """
    print('Generating gmsh input file...\n')
    #formalities and error checks
    if geom_input is not None: 
        if not isinstance(geom_input,dict):
            raise TypeError ("'geom_input' is not a dictionary type object. Dict type is expected for the first argument of genGeoFile")
    else:
        geom_input = {}
    if len(electrodes[0])!=len(electrodes[1]):
        raise ValueError('The length of the electrode x and z arrays does not match')
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
            if key == 'buried': bur_idx.append(i); bu_flag = True
            if key == 'borehole1': bh_flag = True
            if key == 'remote':rem_idx.append(i)
        
        if len(surface_idx)>0:# then surface electrodes are present
            elec_x=np.array(electrodes[0])[surface_idx]
            elec_z=np.array(electrodes[1])[surface_idx]
        elif len(surface_idx)==0 and 'surface' not in geom_input:
            #fail safe if no surface electrodes are present to generate surface topography 
            if not isinstance(geom_input, dict):
                geom_input = {}
            max_idx = np.argmax(electrodes[0])
            min_idx = np.argmin(electrodes[0])
            topo_x = [electrodes[0][min_idx],electrodes[0][max_idx]]
            y_idx = np.argmax(electrodes[1])
            topo_z = [electrodes[1][y_idx]+1,electrodes[1][y_idx]+1] # add one so we don't cut into the buried in electrodes with the mesh
            geom_input['surface'] = [topo_x,topo_z]
            elec_x = np.array([])
            elec_z = np.array([])
            if dp_len == -1 : 
                dp_len = max(electrodes[0]) - min(electrodes[0])
            if doi == -1: 
                doi = min(electrodes[1]) - dp_len
        elif len(surface_idx)==0 and 'surface' in geom_input:
            elec_x = np.array([])
            elec_z = np.array([])
    
    else:
        elec_x = np.array(electrodes[0])
        elec_z = np.array(electrodes[1])

    if 'surface' not in geom_input: # make extra topography points if there are none 
        min_idx = np.argmin(elec_x)
        max_idx = np.argmax(elec_x)   
        topo_x = [elec_x[min_idx] - 5*np.mean(np.diff(np.unique(elec_x))),
                  elec_x[max_idx] + 5*np.mean(np.diff(np.unique(elec_x)))]
        topo_z = [elec_z[min_idx],elec_z[max_idx]]
    else:
        topo_x = geom_input['surface'][0]
        topo_z = geom_input['surface'][1]

    #catch where x coordinates dont change
    if bh_flag or bu_flag:
        elecspacing = (min(electrodes[0]) - max(electrodes[0]))/len(electrodes)
        if max(topo_x)-min(topo_x) < 0.2*elecspacing and len(topo_x) != 1: # they have the same x coordinate 
            topo_x = [min(electrodes[0])-5,max(electrodes[0])+5]
                  
    if dp_len == -1 and len(elec_x)>0:#compute maximum dipole length
        dp_len = abs(np.max(elec_x) - np.min(elec_x))
    elif dp_len == 0 or dp_len==-1:#there is no change dipole length, maybe due to 1 x coordinate 
        dp_len = abs(np.max(electrodes[1]) - np.min(electrodes[1])) # use z coordinates, 
    if dp_len == 0: #if its still 0 
        dp_len = 5 # insert obligitory value here
    if doi == -1:#then set to a default 
        doi = np.max(elec_z) - (np.min(elec_z) - abs(np.max(elec_x) - np.min(elec_x))/2)
    if doi >= np.min(electrodes[1]):
        warnings.warn('DOI is shallower than minimum z coordinate, therefore lowering the DOI by 5 units below lowest electrode')
        doi= np.min(electrodes[1]) - 5 
    print("doi in gmshWrap.py: %f"%doi)
    print("dp_len in gmshWrap.py: %f"%dp_len)
    
    if cl == -1:
        if bh_flag or bu_flag:
            cl = abs(np.mean(np.diff(electrodes[1]))/2)
        else:
            cl = abs(np.mean(np.diff(elec_x))/2)
            
    if len(topo_x) != len(topo_z):
        raise ValueError("topography x and y arrays are not the same length!")
    if len(elec_x) != len(elec_z):
        raise ValueError("electrode x and y arrays are not the same length!")
    
    if file_path.find('.geo')==-1:
        file_path=file_path+'.geo'#add file extension if not specified already
    
    #start to write the file  
    fh = open(file_path,'w') #file handle
    
    fh.write("//2D mesh script for ResIPy (run the following in gmsh to generate a triangular mesh with topograpghy)\n")
    fh.write("cl=%.2f;//define characteristic length\n" %cl)
    fh.write("//Define surface points\n")
    #we have surface topograpghy, and electrode positions to make use of here:
    x_pts=np.append(topo_x,elec_x)#append our electrode positions to our topograpghy 
    y_pts=np.append(topo_z,elec_z)
    flag=['topography point']*len(topo_x)
    flag=flag+(['electrode']*len(elec_x))   
    
    # increase the size of the fine mesh region on both side of the survey
#    e_left = np.min(elec_x) - 4*np.median(np.abs(np.diff(elec_x)))
#    e_right = np.max(elec_x) + 4
    
    #deal with end case electrodes, check max topo points are outside survey bounds 
    try:
        min_idx = np.argmin(elec_x)
        max_idx = np.argmax(elec_x)
        if min(elec_x) == min(x_pts):
            x_pts = np.append(x_pts,elec_x[min_idx] - 5*np.mean(np.abs(np.diff(elec_x)))) # in this case extend the survey bounds beyond the first electrode 
            y_pts = np.append(y_pts,elec_z[min_idx])
            flag.append('topography point')#add a flag
            
        if max(elec_x) == max(x_pts):
            x_pts = np.append(x_pts,elec_x[max_idx] + 5*np.mean(np.abs(np.diff(elec_x))))
            y_pts = np.append(y_pts,elec_z[max_idx])
            flag.append('topography point')
    
    except Exception as e: # then there are no surface electrodes, in which case
        print('Error in assigning additional points on each side of transect: ', e)
        min_idx = np.argmin(electrodes[0])
        max_idx = np.argmax(electrodes[0])
        if min(electrodes[0]) == min(x_pts):
            x_pts = np.append(x_pts,electrodes[0][min_idx] - 5*np.mean(np.diff(electrodes[0]))) # in this case extend the survey bounds beyond the first electrode 
            y_pts = np.append(y_pts,max(electrodes[1])+1)
            flag.append('topography point')#add a flag
            
        if max(electrodes[0]) == max(x_pts):
            x_pts = np.append(x_pts,electrodes[0][max_idx] + 5*np.mean(np.diff(electrodes[0])))
            y_pts = np.append(y_pts,max(electrodes[1])+1)
            flag.append('topography point')
            
    #catch an error where the fine mesh region will truncate where the electrodes are        
    if min(elec_x) < min(x_pts):
        raise ValueError("The minimum X coordinate value for the surface of the mesh must be smaller than the minimum electrode X coordinate")
    if max(elec_x) > max(x_pts):
        raise ValueError("The maximum X coordinate value for the surface of the mesh must be greater than the maximum electrode X coordinate") 
    
    idx=np.argsort(x_pts) # now resort the x positions in ascending order
    x_pts=x_pts[idx] # compute sorted arrays
    y_pts=y_pts[idx]
    z_pts=np.zeros((len(x_pts),1))#for now let z = 0
    #add flags which distinguish what each point is 
    flag_sort=[flag[i] for i in idx]
    
    elec_x_cache = x_pts[np.array(flag_sort)=='electrode']
    elec_z_cache = y_pts[np.array(flag_sort)=='electrode']
    
    #we need protection against repeated points, as this will throw up an error in R2 when it comes calculating element areas
    cache_idx=[]
    for i in range(len(x_pts)-1):
        if x_pts[i]==x_pts[i+1] and y_pts[i]==y_pts[i+1]:
            cache_idx.append(i)
     
    if len(cache_idx)>0:
        warnings.warn("Duplicated surface and electrode positions were detected, R2 inversion likley to fail due to the presence of elements with 0 area.")
        #if duplicated points were dectected we should remove them?? - now disabled 
#        x_pts=np.delete(x_pts,cache_idx)#deletes first instance of the duplicate       
#        y_pts=np.delete(y_pts,cache_idx)
#        print("%i duplicated coordinate(s) were deleted" %len(cache_idx))
#        for k in range(len(cache_idx)):#what are we actually deleting?
#            if flag_sort[cache_idx[k]] == 'electrode':
#                flag_sort[cache_idx[k]+1] = 'electrode'
#                #this overwrites the flag string to say that this point is an electrode
#                #otherwise we'll get a mismatch between the number of electrodes and mesh nodes assigned to the electrodes
#        flag_sort=np.delete(flag_sort,cache_idx).tolist()
    
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
    
    #reflect surface topography at base of fine mesh area, should be a smoothed version of the topography
    x_base = x_pts
#    z_base = moving_average(y_pts - abs(doi),N=5) # compute the depth to the points at the base of the survey
    z_base = y_pts - abs(doi)
    
    # check that the depth lowermost point of the fine mesh region doesnt intersect the surface
    if np.max(z_base) > np.min(electrodes[1]):
        warnings.warn("The depth of investigation is above the the minium z coordinate of electrodes, mesh likely to be buggy!", Warning)   
    
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
#    fh.write("Plane Surface(1) = {1};//Fine mesh region surface\n")
    
    #now extend boundaries beyond flanks of survey area (so generate your Neummon boundary)
    fh.write("\n//Background region (Neumann boundary) points\n")
    cl_factor2=25*cl_factor#characteristic length multipleier for Neumann boundary 
    cln=cl*cl_factor2#assign new cl, this is so mesh elements get larger from the main model
    fh.write("cln=%.2f;//characteristic length for background region\n" %cln)
    #Background region propeties, follow rule of thumb that background should be 5*largest dipole 
    flank=5*dp_len
    b_max_depth=-abs(doi)-(3*dp_len)#background max depth
    #add Neumann boundaries on left hand side
    n_pnt_cache=[0,0,0,0]#cache for the indexes of the neumon boundary points 
    tot_pnts=tot_pnts+1;n_pnt_cache[0]=tot_pnts
    fh.write("Point(%i) = {%.2f,%.2f,%.2f,cln};//far left upper point\n"%(tot_pnts,x_pts[0]-flank,y_pts[0],z_pts[0]))
    tot_pnts=tot_pnts+1;n_pnt_cache[1]=tot_pnts
    fh.write("Point(%i) = {%.2f,%.2f,%.2f,cln};//far left lower point\n"%(tot_pnts,x_pts[0]-flank,b_max_depth,z_pts[0]))
    #add Neumann boundary points on right hand side
    tot_pnts=tot_pnts+1;n_pnt_cache[2]=tot_pnts
    fh.write("Point(%i) = {%.2f,%.2f,%.2f,cln};//far right upper point\n"%(tot_pnts,x_pts[-1]+flank,y_pts[-1],z_pts[-1]))
    tot_pnts=tot_pnts+1;n_pnt_cache[3]=tot_pnts
    fh.write("Point(%i) = {%.2f,%.2f,%.2f,cln};//far right lower point\n"%(tot_pnts,x_pts[-1]+flank,b_max_depth,z_pts[-1]))
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
    coarse_msh_loop = n_ln_cache + [end_ln_cache[1]] + basal_ln_cache_flipped + [-1*end_ln_cache[0]]
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
    #print('probing for boundaries and other additions to the mesh')
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
                
            bhx = electrodes[0][bh_idx]#get borehole coordinate information
            bhy = electrodes[1][bh_idx]
            #cache the x and y coordinates 
            elec_x_cache = np.append(elec_x_cache,bhx)
            elec_z_cache = np.append(elec_z_cache,bhy)
            e_pt_idx = [0] *len(bhx)
            fh.write("// string electrodes for borehole %i\n"%count)
            for k in range(len(bhx)):
                no_pts += 1 
                fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl};//borehole %i electrode\n"%(no_pts,bhx[k],bhy[k],0,count))
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
    
    no_plane = 1 # number of plane surfaces so far (actually two)
    fh.write("\n//Adding polygons\n")
    line_loops = []
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
            fh.write("//make lines forming polygons from line loops\n")
            no_lin += 1
#            fh.write("Line{%s} In Surface{1};\n"%str(line_idx).strip('[').strip(']'))
            fh.write("Line Loop(%i) = {%s};\n"%(no_lin,str(line_idx).strip('[').strip(']')))
            line_loops.append(no_lin)
            no_plane +=1
            fh.write("Plane Surface(%i) = {%i};\n"%(no_plane,no_lin))
            fh.write("Physical Surface (%i) = {%i};\n"%(no_plane,no_plane))
            
        except KeyError:
            fh.write("//end of polygons.\n")
            print('%i polygons added to input file'%(count-1))
            break  
        
    line_loops = np.array(line_loops)
    no_plane += 1
    if len(line_loops) > 0:
        fh.write("Plane Surface(%i) = {%s, 1};//Fine mesh region surface\n"%(no_plane, ', '.join(line_loops.astype(str))))
    else:
        fh.write("Plane Surface(%i) = {1};//Fine mesh region surface\n"%(no_plane))
    fh.write("\n//Make a physical surface\n")
    fh.write("Physical Surface(1) = {%i, 1};\n"%no_plane)
    
    
    #add buried electrodes? (added after as we need Surface 1 to be defined)       
    if bu_flag:
        fh.write("\n//Buried electrodes \n")  
        buried_x = np.array(electrodes[0])[bur_idx]#get buried electrode coordinate information
        buried_z = np.array(electrodes[1])[bur_idx]
        buried_y = [0]*len(buried_x)
        e_pt_idx = [0]*len(buried_x)
        for k in range(len(buried_x)):
            no_pts += 1 
            fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl};//buried electrode %i\n"%(no_pts,buried_x[k],buried_z[k],buried_y[k],k+1))
            e_pt_idx[k] = no_pts
        
        node_pos = np.append(node_pos,e_pt_idx) #add borehole nodes to electrode node positions 
        fh.write("Point{%s} In Surface{%i};\n"%(str(e_pt_idx).strip('[').strip(']'), no_plane))
        fh.write("//end of buried electrodes.\n")
        elec_x_cache = np.append(elec_x_cache,buried_x)
        elec_z_cache = np.append(elec_z_cache,buried_z)
        
    if len(rem_idx)>0: #then we need to add remote electrodes to the mesh
        fh.write("\n//Remote electrodes \n")
        e_pt_idx = [0]*len(rem_idx)
        #put points in lower left hand corner of the outer region
        #if an electrode does have the remote flag, then chances are it has 
        #a bogus coordinate associated with it (ie. -99999)
        for k in range(len(rem_idx)):
            no_pts += 1
            remote_x = x_pts[0]-flank + 10*np.random.rand()
            remote_z = b_max_depth  + 10*np.random.rand()
            fh.write("Point(%i) = {%.2f,%.2f,%.2f,cl2};//remote electrode\n"%(no_pts,remote_x,remote_z,z_pts[0]))
            e_pt_idx[k] = no_pts
            elec_x_cache = np.append(elec_x_cache,electrodes[0][rem_idx[k]])
            elec_z_cache = np.append(elec_z_cache,electrodes[1][rem_idx[k]])
            
        node_pos = np.append(node_pos,e_pt_idx) #add remote electrode nodes to electrode node positions 
        fh.write("Point{%s} In Surface{1};\n"%(str(e_pt_idx).strip('[').strip(']')))
        fh.write('//End of remote electrodes.\n')
        
    
    fh.write("\n//Adding boundaries\n")
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
    
    if len(node_pos) != len(elec_x_cache):
        warnings.warn("looks like something has gone wrong with node orderings, total x != total nodes.")
    
    #sort node ordering back into the original input ordering    
    original_x = np.array(electrodes[0])
    original_z = np.array(electrodes[1])
    #find the original indexes  
    original_idx = [0]*len(node_pos)
    for i in range(len(node_pos)):
        idx = (elec_x_cache == original_x[i]) & (elec_z_cache == original_z[i])
        idx = idx.tolist()
        original_idx[i] = idx.index(True)

    ordered_node_pos = node_pos[original_idx].astype(int)
    
    return ordered_node_pos 

#%% parse a .msh file
def msh_parse(file_path):
    """
    Converts a 2d gmsh mesh file into a mesh class used in pyR2
    
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
    print("ignoring %i non-triangle elements in the mesh file, as they are not required for R2"%ignored_elements)
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
    try:
        if min(areas)==0:
            warnings.warn("elements with no area have been detected in 'mesh.dat', inversion with R2 unlikey to work!" )
    except ValueError:#if mesh hasnt been read in this is where the error occurs 
        raise Exception("It looks like no elements have read into pyR2, its likley gmsh has failed to produced a stable mesh. Consider checking the mesh input (.geo) file.")
            
    print("%i element node orderings had to be corrected because they were found to be orientated clockwise\n"%num_corrected)
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
            'node_y':z_coord,#y coordinates of nodes - nb swapped around as meshTools sees z as elevation
            'node_z':y_coord,#z coordinates of nodes 
            'node_id':node_num,#node id number 
            'elm_id':np.arange(1,real_no_elements+1,1),#element id number 
            'num_elm_nodes':3,#number of points which make an element
            'node_data':node_dump,#nodes of element vertices
            'elm_centre':(centriod_x,[0]*len(centriod_x),centriod_y),#centre of elements (x,y)
            'elm_area':areas,
            'cell_type':[5],
            'parameters':phys_entity,#the values of the attributes given to each cell 
            'parameter_title':'material',
            'dict_type':'mesh_info',
            'original_file_path':file_path} 
    
    #we could return a mesh object here, but a dictionary is easier to debug with spyder, 
    #also we'd need to import the mesh class, and its not a good idea to have modules
    #importing each other, as meshTools has a dependency on gmshWrap. 
        
#%% 2D whole space 

def gen_2d_whole_space(electrodes, padding = 20, electrode_type = None, geom_input = None,
                       file_path='mesh.geo',cl=-1,cl_factor=50,doi=None,dp_len=None):
    """
    writes a gmsh .geo for a 2D whole space. Ignores the type of electrode. 
    
    Parameters
    ----------
    electrodes: array like
        first column/list is the x coordinates of electrode positions, second column
        is the elevation
    padding: float, optional
        Padding in percent on the size the fine mesh region extent. Must be bigger than 0.
        
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
    doi: N/A
        Parameter is not used. Needed for keyword argument compatiblity with 
        genGeoFile.
    dp_len: N/A
        Parameter is not used. Needed for keyword argument compatiblity with 
        genGeoFile.            
    
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

    #### TODO: search through each set of points and check for repeats ?
    """
    
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
        
    def find_dist(elec_x,elec_y,elec_z): # find maximum and minimum electrode spacings 
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
    
    if cl==-1:
        dist_sort = np.unique(find_dist(elec_x,[0]*len(elec_x),elec_z))
        cl = dist_sort[0]/2 # characteristic length is 1/2 the minimum electrode distance
        
    fh = open(file_path,'w') #file handle
    
    fh.write("//Gmsh wrapper code version 1.0 (run the following in gmsh to generate a triangular mesh for 2D whole space)\n")
    fh.write("//2D mesh coordinates\n")
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
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl_pad};\n"%(no_pts, max_x, max_z, 0))
    no_pts += 1
    loop_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl_pad};\n"%(no_pts, max_x, min_z, 0))
    no_pts += 1
    loop_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl_pad};\n"%(no_pts, min_x, min_z, 0))
    no_pts += 1
    loop_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl_pad};\n"%(no_pts, min_x, max_z, 0))
    
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
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl2};\n"%(no_pts, max_x+flank_x, max_z+flank_z, 0))
    no_pts += 1
    nmn_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl2};\n"%(no_pts, max_x+flank_x, min_z-flank_z, 0))
    no_pts += 1
    nmn_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl2};\n"%(no_pts, min_x-flank_x, min_z-flank_z, 0))
    no_pts += 1
    nmn_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl2};\n"%(no_pts, min_x-flank_x, max_z+flank_z, 0))
    
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
    node_pos = [0]*len(elec_x)
    for i in range(len(elec_x)):
        no_pts += 1
        node_pos[i] = no_pts
        fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl};\n"%(no_pts, elec_x[i], elec_z[i], 0))
        fh.write("Point{%i} In Surface{1};\n"%(no_pts))# put the point surface
    fh.write("//End electrodes\n")
    
    fh.write("\n//Adding polygons?\n")
    no_lin=no_lns
    no_plane=0
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
    
    fh.close()
    print("writing .geo to file completed, save location:\n%s\n"%os.getcwd())
    return np.array(node_pos)

#%% 3D half space 

def box_3d(electrodes, padding=20, doi=-1, file_path='mesh3d.geo',
           cl=-1, cl_factor=8, cln_factor=100, dp_len=-1, mesh_refinement=None):
    """
    writes a gmsh .geo for a 3D half space with no topography. Ignores the type of electrode. 
    Z coordinates should be given as depth below the surface! If Z != 0 then its assumed that the
    electrode is buried. 
    
    Parameters
    ----------
    electrodes: list of array likes
        first column/list is the x coordinates of electrode positions, second column
        is the elevation. Z coordinates must normalised to the flat surface if given.ie. 
        z is the depth below the surface. 
    padding: float, optional
        Padding in percent on the size the fine mesh region extent. Must be bigger than 0.
    doi: float, optional 
        Depth of investigation of the survey. 
    file_path: string, optional 
        name of the generated gmsh file (can include file path also) (optional)
    cl: float, optional
        characteristic length (optional), essentially describes how big the nodes 
        assocaited elements will be. Usually no bigger than 5. If set as -1 (default)
        a characteristic length 1/4 the minimum electrode spacing is computed.
    cl_factor: float, optional 
        This allows for tuning of the incremental size increase with depth in the 
        mesh, usually set to 2 such that the elements at the DOI are twice as big as those
        at the surface. The reasoning for this is because the sensitivity of ERT drops
        off with depth. 
    cln_factor: float, optional
        Factor applied to the characteristic length for fine mesh region to compute
        a characteristic length for background (nuemmon) region
    mesh_refinement: list of array likes 
        Coordinates for discrete points in the mesh. 
    
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
    
    def find_dist(elec_x,elec_y,elec_z): # find maximum and minimum electrode spacings 
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
    
    triggered=False
    if cl==-1:
        dist_sort = np.unique(find_dist(elec_x,elec_y,elec_z))
        cl = dist_sort[0]/2 # characteristic length is 1/2 the minimum electrode distance
        triggered = True
        #print('cl = %f'%cl)  
        
    if dp_len == -1: # compute largest possible dipole length 
        if not triggered:# Avoid recalculating if done already 
            dist_sort = np.unique(find_dist(elec_x,elec_y,elec_z))
        dp_len = dist_sort[-1] # maximum possible dipole length
    
    if doi == -1: # compute depth of investigation
        try:
            doi = dist_sort[-1]/3 # maximum possible dipole length / 3
        except NameError: # putting in try to avoid recalculating if done already 
            dist_sort = np.unique(find_dist(elec_x,elec_y,elec_z))
            doi = dist_sort[-1]/3
    if doi >= np.min(elec_z):
        warnings.warn('depth of investigation is shaller than lowest electrode, adding 5 units below lowest electrode')
        doi = np.min(elec_z) - 5 
        
    print('doi in gmshWrap.py: %f'%doi)
    ### start to write to file ### 
    fh = open(file_path,'w') #file handle
    
    fh.write("//3D half space problem mesh for ResIPy - no topography\n")
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
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl};\n"%(no_pts, max_x, max_y, 0))
    no_pts += 1
    loop_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl};\n"%(no_pts, max_x, min_y, 0))
    no_pts += 1
    loop_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl};\n"%(no_pts, min_x, min_y, 0))
    no_pts += 1
    loop_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl};\n"%(no_pts, min_x, max_y, 0))
    
    #add line loop
    no_lns = 0 
    for i in range(4):
        no_lns += 1 
        if i == 3:
            fh.write("Line(%i) = {%i,%i};\n"%(no_lns,loop_pt_idx[i],loop_pt_idx[0]))
        else:
            fh.write("Line(%i) = {%i,%i};\n"%(no_lns,loop_pt_idx[i],loop_pt_idx[i+1]))
    
    doi = abs(doi)
    #add points below surface to make a rectangular volume  
    fh.write("cl2=cl*%.2f;//define characteristic length for base of fine mesh region\n" %cl_factor)       
    no_pts += 1
    loop2_pt_idx=[no_pts]
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl2};\n"%(no_pts, max_x, max_y, -doi))
    no_pts += 1
    loop2_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl2};\n"%(no_pts, max_x, min_y, -doi))
    no_pts += 1
    loop2_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl2};\n"%(no_pts, min_x, min_y, -doi))
    no_pts += 1
    loop2_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl2};\n"%(no_pts, min_x, max_y, -doi))
    
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
    flank_z = 3*abs(doi)
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
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cln};\n"%(no_pts, max_x+flank_x, max_y+flank_y, -doi - flank_z))
    no_pts += 1
    nmn2_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cln};\n"%(no_pts, max_x+flank_x, min_y-flank_y, -doi - flank_z))
    no_pts += 1
    nmn2_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cln};\n"%(no_pts, min_x-flank_x, min_y-flank_y, -doi - flank_z))
    no_pts += 1
    nmn2_pt_idx.append(no_pts)
    fh.write("Point (%i) = {%.2f,%.2f,%.2f, cln};\n"%(no_pts, min_x-flank_x, max_y+flank_y, -doi - flank_z))
    
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
    
    fh.write("//Electrode positions.\n")
    node_pos = [0]*len(elec_x)
    for i in range(len(elec_x)):
        no_pts += 1
        node_pos[i] = no_pts
        if elec_z[i] == 0:
            fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl};\n"%(no_pts, elec_x[i], elec_y[i], 0))
            fh.write("Point{%i} In Surface{1};\n"%(no_pts))# put the point surface
        else:
            fh.write("Point (%i) = {%.2f,%.2f,%.2f, cl};\n"%(no_pts, elec_x[i], elec_y[i], elec_z[i]))
            fh.write("Point{%i} In Volume{1};\n"%(no_pts))# put the point in volume 
    fh.write("//End electrodes\n")
    print("writing .geo to file completed, save location:\n%s\n"%os.getcwd())
    return np.array(node_pos) 
    
#%% parse a .msh file
def msh_parse_3d(file_path):
    """
    Converts a 3d gmsh mesh file into a mesh class used in pyR2
    
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
    count = 0 
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
            no_elements=int(dump[i+1])#number of elements
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
    node3 = []#last node of triangle
    node4 = [] 
    npere = 4 
    ignored_elements=0#count the number of ignored elements
    for i in range(element_start,element_end):
        line_data=[int(k) for k in dump[i].split()]
        if line_data[1]==4:# then its the right element type!
        #convert string info into floats and cache data
            nat_elm_num.append(line_data[0])
            elm_type.append(line_data[1]) 
            number_of_tags.append(line_data[2]) 
            phys_entity.append(line_data[3]) 
            elem_entity.append(line_data[4]) 
            node1.append(line_data[5]-1) 
            node2.append(line_data[6]-1) 
            node3.append(line_data[7]-1)
            node4.append(line_data[8]-1)
        else:
            ignored_elements += 1
    print("ignoring %i non-tetrahedra elements in the mesh file, as they are not required for R3t"%ignored_elements)
    real_no_elements=len(nat_elm_num) #'real' number of elements that we actaully want
    
    #compute element centres 
    centriod_x=[]
    centriod_y=[]
    centriod_z=[]
    areas=[]
    for i in range(real_no_elements):
        n1=(x_coord[node1[i]],y_coord[node1[i]],z_coord[node1[i]])#define node coordinates
        n2=(x_coord[node2[i]],y_coord[node2[i]],z_coord[node2[i]])#we have to take 1 off here cos of how python indexes lists and tuples
        n3=(x_coord[node3[i]],y_coord[node3[i]],z_coord[node3[i]])
        n4=(x_coord[node4[i]],y_coord[node4[i]],z_coord[node4[i]])
        centriod_x.append(sum((n1[0],n2[0],n3[0],n4[0]))/npere)
        centriod_y.append(sum((n1[1],n2[1],n3[1],n4[1]))/npere)
        centriod_z.append(sum((n1[2],n2[2],n3[2],n4[2]))/npere)
        
    node_dump = [node1,node2,node3,node4]
    
    mesh_dict = {'num_elms':real_no_elements,
            'num_nodes':no_nodes,
    #        'num_regions':0,
    #        'element_ranges':assctns,
            'dump':dump,      
            'node_x':x_coord,#x coordinates of nodes 
            'node_y':y_coord,#y coordinates of nodes
            'node_z':z_coord,#z coordinates of nodes 
            'node_id':node_num,#node id number 
            'elm_id':np.arange(1,real_no_elements+1,1),#element id number 
            'num_elm_nodes':3,#number of points which make an element
            'node_data':node_dump,#nodes of element vertices
            'elm_centre':[centriod_x,centriod_y,centriod_z],#centre of elements (x,y,z)
            'elm_area':areas,
            'cell_type':[10],
            'parameters':phys_entity,#the values of the attributes given to each cell 
            'parameter_title':'material',
            'dict_type':'mesh_info',
            'original_file_path':file_path} 
    return mesh_dict 
