# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 08:19:06 2018
Import a vtk file with an unstructured grid 
@author: jamyd91
functions: 
    tri_cent
    vtk_import
"""
#import tkinter as tk
#from tkinter import filedialog

#%% triangle centriod 
def tri_cent(p,q,r):#code expects points as p=(x,y) and so on ... (counter clockwise prefered)
    Xm=(p[0]+q[0])/2
    Ym=(p[1]+q[1])/2
    k=2/3
    Xc=r[0]+(k*(Xm-r[0]))
    Yc=r[1]+(k*(Ym-r[1]))
    return(Xc,Yc)
    
#%% import a vtk file 
def vtk_import(file_path='ask_to_open',parameter_title='default'):
    #Converts a gmsh mesh file into a mesh.dat file needed for R2. 
#INPUT:
    #file_path - file path to mesh file. note that a error will occur if the file format is not as expected
    #save_path - leave this as default to save the file in the working directory, make this 'ask_to_open' to open a dialogue box, else enter a custom file path.
    #parameter_title - name of the parameter table in the vtk file, if left as default the first look up table found will be returned 
#OUTPUT: 
    #dictionary with some 'useful' info about the mesh
###############################################################################
#    if file_path=='ask_to_open':#use a dialogue box to open a file
#        print("please select the vtk file to import using the pop up dialogue box. \n")
#        root=tk.Tk()
#        root.withdraw()
#        file_path=filedialog.askopenfilename(title='Select mesh file',filetypes=(("VTK files","*.vtk"),("all files","*.*")))#
    #open the selected file for reading
    fid=open(file_path,'r')
    print("importing vtk (2D mesh) file into python workspace...")
    #read in header info and perform checks to make sure things are as expected
    vtk_ver=fid.readline().strip()#read first line
    if vtk_ver.find('vtk')==-1:
        raise ImportError("Unexpected file type... ")
    elif vtk_ver.find('3.0')==-1:#not the development version for this code
        print("Warning: vtk manipulation code was developed for vtk datafile version 3.0, unexpected behaviour may occur")
    title=fid.readline().strip()#read line 2
    format_type=fid.readline().strip()#read line 3
    if format_type=='BINARY':
        raise ImportError("expected ASCII type file format, not binary")
    dataset_type=fid.readline().strip().split()#read line 4
    if dataset_type[1]!='UNSTRUCTURED_GRID':
        print("Warning: code intended to deal with an 'UNSTRUCTURED_GRID' data type not %s"%dataset_type[1])
    #read node data
    print("importing mesh nodes...")
    node_info=fid.readline().strip().split()#read line 5
    no_nodes=int(node_info[1])
    #now read in node data
    x_coord=[]#make lists for each of the relevant parameters for each node
    y_coord=[]
    z_coord=[]
    node_num=[]
    for i in range(no_nodes):
        coord_data=fid.readline().strip().split()
        x_coord.append(float(coord_data[0]))
        y_coord.append(float(coord_data[1]))
        z_coord.append(float(coord_data[2]))
        node_num.append(i)
    #now read in element data
    print("importing mesh element info...")
    elm_info=fid.readline().strip().split()#read line with cell data
    no_elms=int(elm_info[1])
    no_pts=[]
    node1=[]
    node2=[]
    node3=[]
    elm_num=[]
    centriod_x=[]
    centriod_y=[]
    areas=[]
    #import element data ... expects triangles
    for i in range(no_elms):
        elm_data=fid.readline().strip().split()
        if int(elm_data[0])!=3:
            raise ImportError("non-triangle elements detected! Aborting vtk import.")
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
    centriod=(centriod_x,centriod_y)#centres of each element in form (x...,y...)
    #now for final part of file - cell type info
    cell_type_data=fid.readline().strip()
    cell_type=fid.readline().strip().split()
    _=fid.readline()#read point data line
    _=fid.readline()#read cell data line ... i'm not sure why these need to be repeated, must be for the table lookup process
    cell_attributes=fid.readlines()#reads the last portion of the file
    #finished reading the file
    fid.close()
    print("reading cell attributes...")
    # read through cell attributes to find the relevant parameter table?
    if parameter_title=='default' and title=='Output from R2':    
        parameter_title='Resistivity(Ohm-m)'# the name of title if the output is from R2
        do_find=1
    elif parameter_title == 'n/a':#dont bother looking for attributes
        do_find=0
    elif parameter_title=='default':
        do_find=2
    else:
        do_find=1
    #now that conditions for finding a parameter table have been decided... 
    if do_find==1:
        for i in range(len(cell_attributes)):
            probe=cell_attributes[i].split()
            if probe[1]==parameter_title:
               #then the following line should read "LOOKUP_TABLE default"
               check=cell_attributes[i+1]
               print("identified relevant table for element attributes...")
               indx=i+2
               break
            if i==range(len(cell_attributes)):
               print("Warning: could not find relevant table for element attributes! Make sure you havent made a mistake with table name in the VTK file. \n")
               indx=3
        values=[float(k) for k in cell_attributes[indx].split()]
    elif do_find==2:
        if len(cell_attributes)>=3:
            probe=cell_attributes[1].split()
            parameter_title=probe[1]
            values=[float(k) for k in cell_attributes[3].split()]
        else:
            values='n/a'    
    elif do_find==0:
        values='n/a'
#need two options here, either find depth or find if the elements lie in a certain region
    print("finished importing mesh.\n")
#return information in a dictionary: 
    return {'num_nodes':no_nodes,#number of nodes
            'num_elms':no_elms,#number of elements 
            'node_x':x_coord,#x coordinates of nodes 
            'node_y':y_coord,#y coordinates of nodes
            'node_z':z_coord,#z coordinates of nodes 
            'node_id':node_num,#node id number 
            'elm_id':elm_num,#element id number 
            'num_elm_nodes':no_pts,#number of points which make an element
            'node_data':(node1,node2,node3),#nodes of element vertices
            'elm_centre':centriod,#centre of elements (x,y)
            'elm_area':areas,#area of each element
            'cell_type':cell_type,
            'parameters':values,#the values of the attributes given to each cell 
            'parameter_title':parameter_title,
            'cell_attribute_dump':cell_attributes,
            'dict_type':'mesh_info',
            'original_file_path':file_path} 
    