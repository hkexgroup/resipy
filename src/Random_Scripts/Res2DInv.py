#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15th 2019
Handles importing external data files into the R2 api
@author: Pedro Concha

Currently supports: 
    res2dinv input files
"""

import numpy as np
import pandas as pd
import os 


#%% function to compute geometric factor - written by Pedro Concha 
def geom_fac_1(C1,C2,P1,P2):
    Rc1p1 = np.abs(P1 - C1)
    Rc2p1 = np.abs(C2 - P1)
    Rc1p2 = np.abs(P2 - C1)
    Rc2p2 = np.abs(C2 - P2)
    
    denom = (1/Rc1p1) - (1/Rc2p1) - (1/Rc1p2) + (1/Rc2p2)
    k = (2*np.pi)/denom
    return k 

#%% function to read Res2Dinv files - written by Pedro Concha
def res2DinvParser(file_path):
    """
    Reads *.dat files in index data format and general array format.
    The following arrays are supported: Wenner alpha, Beta and Gamma, dipole-dipole,
    Schlumberger.
    
    Parameters
    -----------
    file_path : string 
     string mapping to the res2inv input file 
    
    Returns
    ----------
    elec : np array
        electrode coordinate matrix in the form | x | y | z |
    df: pandas dataframe
        dataframe which holds the electrode numbers for in feild measurements and 
        apparent resistivities (Rho) and transfer resistances 
    
    ## TODO : Pole-Pole, Pole-Dipole and borehole surveys 
    """
    fh = open(file_path,'r')#open file handle for reading
    dump = fh.readlines()#cache file contents into a list
    fh.close()#close file handle, free up resources
    
    #first find loke's General array format information in the file (which is all we need for R2/3)
    fmt_flag = False # format flag 
    err_flag = False # error estimates flag 
    topo_flag = False # topography flag 
    start_0_flag = False
    sur_flag = 0 # survey type
    x_location = 0 # for data points, 0 for first electrode, 1 for mid-point, 2 for surface distance
    ip_flag = 0 # 0 for none IP data, 1 if present
    idx_oi = 0
    line = dump[idx_oi]
    sur_name = line.strip() #name of survey
    idx_oi += 1
    line = dump[idx_oi]
    a_spac = float(line) #electrode spacing
    idx_oi += 1
    line = dump[idx_oi]
    array_type = int(line)
    if array_type in [2,6]:
        raise ImportError("Not supported")      
       
    if array_type in [1,3,4,5,7]: 
        #1:Wenner alfa, 3:DipDip,4:Wenner beta, 5:Wenner gama,
        #7:Schlumberger, 11:General array, 15:Gradient 
        idx_oi = 3
        line = dump[idx_oi]
    elif array_type in (11, 15):
        idx_oi = 6
        line = dump[idx_oi]
    
    num_meas = int(line)
    idx_oi += 1
    line = dump[idx_oi]
    x_location = int(line)
    idx_oi += 1
    line = dump[idx_oi]
    ip_flag = int(line)
    idx_oi += 1
    line = dump[idx_oi]
    
    Pa = []
    x_dump = []
    total_x=np.array(())
    e_idx = []
    total_x = np.append(total_x,0)
    data_dict = {'a':[],'b':[],'m':[],'n':[],'Rho':[],'ip':[],'resist':[],'dev':[]}
    
    for k in range(num_meas):
        line = dump[idx_oi + k]
        vals = line.strip().split()
        if array_type == 1:# Wenner alfa
            if x_location == 0:
                c1 = (float(vals[0]))
                a = (float(vals[1]))
                p1 = c1 + a
                p2 = p1 + a
                c2 = p2 + a
            if x_location == 1:
                mid_point = (float(vals[0]))
                mn = (float(vals[1]))
                p1 = mid_point - mn/2
                c1 = p1 - mn
                p2 = p1 + mn
                c2 = p2 + mn
        elif array_type == 3:#Dipole-Dipole
            if x_location == 0:
                c2 = (float(vals[0]))
                a = (float(vals[1]))
                c1 = c2 + a
                n = (float(vals[2]))
                p1 = c1 + a * n
                p2 = p1 + a
            if x_location == 1:
                mid_point = (float(vals[0]))
                a = (float(vals[1]))
                n = (float(vals[2]))
                c1 = mid_point - n*a/2
                c2 = c1 - a
                p1 = c1 + n*a
                p2 = p1 + a
        elif array_type == 4:# Wenner beta
            if x_location == 0:
                c2 = (float(vals[0]))
                a = (float(vals[1]))
                c1 = c2 + a 
                p1 = c1 + a
                p2 = p1 + a
            if x_location == 1:
                mid_point = (float(vals[0]))
                a = (float(vals[1]))
                c1 = mid_point - a/2
                c2 = c1 - a
                p1 = c1 + a
                p2 = p1 + a
        elif array_type == 5:# Wenner gamma
            if x_location == 0:
                c1 = (float(vals[0]))
                a = (float(vals[1]))
                p1 = c1 + a
                c2 = p1 + a
                p2 = c2 + a
            if x_location == 1:
                mid_point = (float(vals[0]))
                a = (float(vals[1]))
                p1 = mid_point - a/2
                c1 = p1 - a
                c2 = p1 + a
                p2 = c2 + a
        elif array_type == 7:#Schlumberger
            if x_location == 0:
                c1 = (float(vals[0]))
                a = (float(vals[1]))
                n = (float(vals[2]))
                p1 = c1 + n * a
                p2 = p1 + a
                c2 = p2 + n * a
            if x_location == 1:
                mid_point = (float(vals[0]))
                a = (float(vals[1]))
                n = (float(vals[2]))
                p1 = mid_point - a/2
                c1 = p1 - n*a
                p2 = p1 + a
                c2 = p2 + n*a
        elif array_type in (11,15):
           c1 = (float(vals[1]))
           if c1 == 0 or start_0_flag:
               start_0_flag = True
           c2 = (float(vals[3]))
           p1 = (float(vals[5]))
           p2 = (float(vals[7]))
            
        x_dump.append(c1)
        x_dump.append(p1)
        x_dump.append(p2)
        x_dump.append(c2)
        total_x = np.append(total_x,x_dump)
        #convert the x electrode coordinates into indexes?
        ex_pos = np.unique(total_x)
        x_dump.clear()
    
    for k in range(num_meas):
        line = dump[idx_oi + k]
        vals = line.strip().split()            
        if array_type == 1:# Wenner alfa
            Pa.append(float(vals[2]))
            if x_location == 0:
                c1 = (float(vals[0]))
                a = (float(vals[1]))
                p1 = c1 + a
                p2 = p1 + a
                c2 = p2 + a
            if x_location == 1:
                mid_point = (float(vals[0]))
                mn = (float(vals[1]))
                p1 = mid_point - mn/2
                c1 = p1 - mn
                p2 = p1 + mn
                c2 = p2 + mn
            if c1 == 0 or start_0_flag:
                    start_0_flag = True
        elif array_type == 3:#Dipole-Dipole
                Pa.append(float(vals[3]))
                if x_location == 0:
                    c2 = (float(vals[0]))
                    a = (float(vals[1]))
                    c1 = c2 + a
                    n = (float(vals[2]))
                    p1 = c1 + a * n
                    p2 = p1 + a
                if x_location == 1:
                    mid_point = (float(vals[0]))
                    a = (float(vals[1]))
                    n = (float(vals[2]))
                    c1 = mid_point - n*a/2
                    c2 = c1 - a
                    p1 = c1 + n*a
                    p2 = p1 + a
                if c2 == 0 or start_0_flag:
                        start_0_flag = True
        elif array_type == 4:# Wenner beta
                Pa.append(float(vals[2]))
                if x_location == 0:
                    c2 = (float(vals[0]))
                    a = (float(vals[1]))
                    c1 = c2 + a 
                    p1 = c1 + a
                    p2 = p1 + a
                if x_location == 1:
                    mid_point = (float(vals[0]))
                    a = (float(vals[1]))
                    c1 = mid_point - a/2
                    c2 = c1 - a
                    p1 = c1 + a
                    p2 = p1 + a                
                if c2 == 0 or start_0_flag:
                    start_0_flag = True
        elif array_type == 5:# Wenner gamma
                Pa.append(float(vals[2]))
                if x_location == 0:
                    c1 = (float(vals[0]))
                    a = (float(vals[1]))
                    p1 = c1 + a
                    c2 = p1 + a
                    p2 = c2 + a
                if x_location == 1:
                    mid_point = (float(vals[0]))
                    a = (float(vals[1]))
                    p1 = mid_point - a/2
                    c1 = p1 - a
                    c2 = p1 + a
                    p2 = c2 + a
                if c1 == 0 or start_0_flag :
                    start_0_flag = True
        elif array_type == 7:#Schlumberger
                Pa.append(float(vals[3]))
                if x_location == 0:
                    c1 = (float(vals[0]))
                    a = (float(vals[1]))
                    n = (float(vals[2]))
                    p1 = c1 + n * a
                    p2 = p1 + a
                    c2 = p2 + n * a
                if x_location == 1:
                    mid_point = (float(vals[0]))
                    a = (float(vals[1]))
                    n = (float(vals[2]))
                    p1 = mid_point - a/2
                    c1 = p1 - n*a
                    p2 = p1 + a
                    c2 = p2 + n*a
                if c1 == 0 or start_0_flag:
                    start_0_flag = True
        elif array_type in (11,15):
           c1 = (float(vals[1]))
           if c1 == 0 or start_0_flag:
               start_0_flag = True
           c2 = (float(vals[3]))
           p1 = (float(vals[5]))
           p2 = (float(vals[7]))
           Pa.append(float(vals[9]))
            
        x_dump.append(c1)
        x_dump.append(p1)
        x_dump.append(p2)
        x_dump.append(c2)
    
        if start_0_flag:
            e_idx = [np.where(ex_pos == x_dump[i])[0][0] for i in range(4)]
            e_idx= np.add(e_idx, [1, 1, 1, 1])
        else:
            e_idx = [np.where(ex_pos == x_dump[i])[0][0] for i in range(4)]
    
        data_dict['a'].append(e_idx[0])
        data_dict['b'].append(e_idx[3])
        data_dict['m'].append(e_idx[1])
        data_dict['n'].append(e_idx[2])
        #convert apparent resistivity back in to transfer resistance
        if array_type in [1, 7, 11]:
            K = geom_fac_1(c1, c2, p1, p2)
        elif array_type == 3:
            K = np.pi * n*(n + 1)*(n + 2)*a
        elif array_type == 4:
            K = 6 * np.pi * a
        elif array_type == 5:
            K = 3 * np.pi * a
        
        Pt = Pa[k]/K # transfer resistance
        #add apparent and transfer resistances to dictionary
        data_dict['Rho'].append(Pa[k])
        data_dict['resist'].append(Pt)
        data_dict['ip'].append(0)
        if err_flag:
            err_Pa = float(vals[6])
            err_Pt = err_Pa/K
            data_dict['dev'].append(abs(err_Pt))
        else:
            data_dict['dev'].append(0)
    
        x_dump.clear()
                  
    num_elec = len(ex_pos)
    fmt_flag = True
        
    topo_flag_idx = idx_oi + num_meas
    
    if int(dump[topo_flag_idx]) == 2 :#if we have topography then we should read it into the API
        #print("topography flag activated")
        topo_flag = True
        num_elec_topo =  int(dump[topo_flag_idx+1])
        ex_pos_topo=[0]*num_elec_topo
        ey_pos_topo=[0]*num_elec_topo
        ez_pos_topo=[0]*num_elec_topo # actaully we can't have a z coordinate for 2d data so these will remain as zero
        ey_pos=[0]*num_elec
        ez_pos=[0]*num_elec
        for i in range(num_elec_topo):
            ex_pos_topo[i] = float(dump[topo_flag_idx+2+i].strip().split()[0])
            ey_pos_topo[i] = float(dump[topo_flag_idx+2+i].strip().split()[1])
            
        for i in range(num_elec):
            for j in range(num_elec_topo):
                if ex_pos[i] == ex_pos_topo[j]:
                    ey_pos[i] = ey_pos_topo[j]
        #print(ex_pos,ey_pos)
        elec = np.column_stack((ex_pos,ey_pos,ez_pos))
              
       
    #add some protection against a dodgey file 
    if fmt_flag is False:
        raise ImportError("Error importing res2dinv input file:"+file_path+"\n the file is either unrecognised or unsupported")        
    
    #now we have indexed electrode coordinates in ex_pos :) 
    if not topo_flag: # then we dont have any topography and the electrode positions are simply given by thier x coordinates
        ey_pos=[0]*num_elec
        ez_pos=[0]*num_elec  
        elec = np.column_stack((ex_pos,ey_pos,ez_pos))
       
    df = pd.DataFrame(data=data_dict) # make a data frame from dictionary
    df = df[['a','b','m','n','Rho','dev','ip','resist']] # reorder columns to be consistent with the syscal parser

    return elec,df

# Test of res2Dinv files
    #%% ########## Load ER data ################
# Define directory where example scripts are located 
work_dir = os.path.abspath(r'/Volumes/GoogleDrive/My Drive/Rifle/Data/9_10/Tracer_Injection_tracking/res2dinv')
#res2Dinv_file = 'LANDFILL.DAT'# Wenner alfa example
#res2Dinv_file = 'BLOCKDIP.DAT'# Dipole_Dipole
#res2Dinv_file = 'beta.dat'# Wenner Beta
#res2Dinv_file = 'GAMMA.DAT'# Wenner Gamma
res2Dinv_file = '1010am.DAT'# Wenner_Schlumberger
#res2Dinv_file = 'MIXED.DAT'# General Array fomat
#res2Dinv_file = 'GLADOE2.DAT'# Topography data for index based format files
file_path = os.path.join(work_dir,res2Dinv_file)

if os.path.exists(file_path):
    elec,df = res2DinvParser(file_path)