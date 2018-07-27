#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  1 11:23:23 2018
Handles importing external data files into the R2 api
@author: jkl 

Currently supports: 
    syscal files
    res2dinv input files 
"""

import numpy as np
import pandas as pd

#%% function to compute geometric factor - Jimmy B 
def geom_fac(C1,C2,P1,P2):
    Rc1p1 = (C1 - P1)
    Rc2p1 = (C2 - P1)
    Rc1p2 = (C1 - P2)
    Rc2p2 = (C2 - P2)
    
    denom = (1/Rc1p1) - (1/Rc2p1) - (1/Rc1p2) + (1/Rc2p2)
    k = (2*np.pi)/denom
    return k 

#%% usual syscal parser
def syscalParser(fname, spacing=None):
        df = pd.read_csv(fname, skipinitialspace=True)
        # delete space at the end and the beginning of columns names
        headers = df.columns
        newheaders = list(map(str.strip, headers)) 
        dico = dict(zip(headers, newheaders))
        df = df.rename(index=str, columns=dico)
        df = df.rename(columns={'Spa.1':'a',
                                'Spa.2':'b',
                                'Spa.3':'m',
                                'Spa.4':'n',
                                'In':'i',
                                'Vp':'vp',
                                'Dev.':'dev',
                                'M':'ip', #M1, M2,...Mn are good for now when importing
                                'Sp':'sp'})
        
        array = df[['a','b','m','n']].values
        if spacing == None:
            spacing = np.unique(np.sort(array.flatten()))[1]
        array = np.round(array/spacing+1).astype(int)
        df[['a','b','m','n']] = array
        df['resist'] = df['vp']/df['i']
        imax = int(np.max(array))
        elec = np.zeros((imax,3))
        elec[:,0] = np.arange(0,imax)*spacing
                
        return elec, df
    
#test code
#elec, df = syscalParser('test/syscalFile.csv')


#%% protocol.dat forward modelling parser

def protocolParser(fname):
    x = np.genfromtxt(fname, skip_header=1)
    df = pd.DataFrame(x, columns=['index','a','b','m','n','appResist','resist'])
    df['ip'] = np.nan
    xElec = np.arange(np.max(df[['a','b','m','n']].values))
    elec = np.zeros((len(xElec),3))
    elec[:,0] = xElec
    return elec, df

# test code
#protocolParser('test/protocolFile.dat')
    
    
#%% parse input for res2inv (.dat file) - Jimmy B. 
#jamyd91@bgs.ac.uk
def res2invInputParser(file_path):
    #biulds an R2 protocal.dat from a res2dinv input . dat file. It looks for the general
    #array format in the .dat file. Note you'll need to do some extra work if you want 
    #to extract the electrode locations as well as the program/schedule information
    
#INPUT:
    #file_path - string mapping to the res2inv input file 
#OUTPUT:
    #string which is the protocal.dat information. 
## TODO : add capacity to recognise errors in the input file (not currently read in)
## TODO : Return electrode coordinates (code already finds the coordinates they just need to be returned)
## TODO : add capacity to read in borehole surveys 
###############################################################################
    fh = open(file_path,'r')#open file handle for reading
    dump = fh.readlines()#cache file contents into a list
    fh.close()#close gile handle, free up RAM
    
    #first find loke's General array format information in the file (which is all we need for R2/3)
    fmt_flag = False # format flag 
    err_flag = False # error estimates flag 
    topo_flag = False # topography flag 
    sur_flag = 0 # survey type --> currently unused . 
    for i,line in enumerate(dump):
        if line.strip() == "General array format":
            #print("found array data")
            idx_oi = i+1 # index of interest
            fmt_flag = True
        #find if errors are present 
        if line.strip() == "Error estimate for data present":
            err_flag = True
            print('errors estimates found in file but are not currently supported')
        
    #add some protection against a dodgey file 
    if fmt_flag is False:
        raise ImportError("Error importing res2dinv input file:"+file_path+"\n the file is either unrecognised or unsupported")        
    
    num_meas = int(dump[3])#number of measurements should be stored on the 4th line of the file
       
    #find topography? 
    topo_flag_idx = 6+num_meas
    if err_flag:#changes if we have an error estimate -_-
        topo_flag_idx = 9+num_meas        
    else:
        topo_flag_idx = 6+num_meas
    
    if int(dump[topo_flag_idx]) == 2 :#if we have topography then we should read it into the API
        #print("topography flag activated")
        topo_flag = True
        num_elec =  int(dump[topo_flag_idx+1])
        ex_pos=[0]*num_elec
        ey_pos=[0]*num_elec
        for i in range(num_elec):
            ex_pos[i] = float(dump[topo_flag_idx+2+i].strip().split(',')[0])
            ey_pos[i] = float(dump[topo_flag_idx+2+i].strip().split(',')[1])
        #print(ex_pos,ey_pos)
              
    #since we dont always have all the electrode indexes we need to determine this
    #idea here is to extract all the x locations from the general array format
    total_x=np.array(())
    #print('finding general array electrode coordinates')
    for k in range(num_meas):
        line = dump[k+idx_oi]
        vals = line.strip().split()
        x_dump = [float(vals[1:5][i].split(',')[0]) for i in range(4)]
        total_x = np.append(total_x,x_dump) # this caches all the x locations 
    #extract the unique x values using numpy
    ex_pos = np.unique(total_x)
    #now we have indexed electrode coordinates in ex_pos :) 
    
    #were having to assume the format of loke's file to be in the forM:
    #no.electrodes | C1 | C2 | P1 | P2 | apparent.resistivity
    #print('computing transfer resistances and reading in electrode indexes')
    protocal = '%i\n'%num_meas # string which holds protocal.dat needed for r2
    for k in range(num_meas):
        line = dump[k+idx_oi]
        vals = line.strip().split()
        x_dump = [float(vals[1:5][i].split(',')[0]) for i in range(4)]
        #convert the x electrode coordinates into indexes?
        e_idx = [np.where(ex_pos == x_dump[i])[0][0] for i in range(4)]
        #convert apparent resistivity back in to transfer resistance
        K = geom_fac(x_dump[0],x_dump[1],x_dump[2],x_dump[3])
        Pa = float(vals[5]) # apparent resistivity value
        Pt = Pa/K # transfer resistance
        protocal += '{}\t{}\t{}\t{}\t'.format(e_idx[0],e_idx[1],e_idx[2],e_idx[3])
        protocal += '{}\n'.format(Pt)
    
    return protocal

