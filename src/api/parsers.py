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
import os 

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
        if 'Spa.1' in headers:
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
        else:
            newheaders = list(map(str.strip, headers)) 
            dico = dict(zip(headers, newheaders))
            df = df.rename(index=str, columns=dico)
            df = df.rename(columns={'xA(m)':'a',
                                    'xB(m)':'b',
                                    'xM(m)':'m',
                                    'xN(m)':'n',
                                    'Dev.':'dev',
                                    'M (mV/V)':'ip',
                                    'SP (mV)':'sp',
                                    'VMN (mV)':'vp',
                                    'IAB (mA)':'i',
                                    'M1 (mV/V)':'M1',
                                    'M2 (mV/V)':'M2',
                                    'M3 (mV/V)':'M3',
                                    'M4 (mV/V)':'M4',
                                    'M5 (mV/V)':'M5',
                                    'M6 (mV/V)':'M6',
                                    'M7 (mV/V)':'M7',
                                    'M8 (mV/V)':'M8',
                                    'M9 (mV/V)':'M9',
                                    'M10 (mV/V)':'M10',
                                    'M11 (mV/V)':'M11',
                                    'M12 (mV/V)':'M12',
                                    'M13 (mV/V)':'M13',
                                    'M14 (mV/V)':'M14',
                                    'M15 (mV/V)':'M15',
                                    'M16 (mV/V)':'M16',
                                    'M17 (mV/V)':'M17',
                                    'M18 (mV/V)':'M18',
                                    'M19 (mV/V)':'M19',
                                    'M20 (mV/V)':'M20',
                                    'TM1 (ms)':'TM1'})
        
        array = df[['a','b','m','n']].values
        arrayMin = np.min(np.unique(np.sort(array.flatten())))
        if arrayMin != 0:
            array -= arrayMin
        espacing = np.unique(np.sort(array.flatten()))[1] - np.unique(np.sort(array.flatten()))[0]
        array = np.round(array/espacing+1).astype(int)
        df[['a','b','m','n']] = array
        df['resist'] = df['vp']/df['i']
        imax = int(np.max(array))
        elec = np.zeros((imax,3))
        if spacing is None:
            spacing = espacing
        elec[:,0] = np.arange(0,imax)*spacing
                
        return elec, df
    
#test code
#elec, df = syscalParser('test/syscalFile.csv')
#elec, df = syscalParser('/media/jkl/data/phd/tmp/projects/ahdb/survey2018-08-14/data/ert/18081401.csv', 0.5)
#print(df[['a','b','m','n']])
#print(elec)

#%% protocol.dat forward modelling parser

def protocolParser(fname):
    x = np.genfromtxt(fname, skip_header=1)
    df = pd.DataFrame(x, columns=['index','a','b','m','n','resist','appResist'])
    df['ip'] = np.nan
    xElec = np.arange(np.max(df[['a','b','m','n']].values))
    elec = np.zeros((len(xElec),3))
    elec[:,0] = xElec
    return elec, df

# test code
#protocolParser('test/protocolFile.dat')

def protocolParserIP(fname): # just for protocol forward output with cR2
    x = np.genfromtxt(fname, skip_header=1)
    df = pd.DataFrame(x, columns=['index','a','b','m','n','resist','ip','appResist'])
    xElec = np.arange(np.max(df[['a','b','m','n']].values))
    elec = np.zeros((len(xElec),3))
    elec[:,0] = xElec
    return elec, df    

#%% PRIME system parser

def primeParser(fname, espacing=None):
    """ Returns data and elec from BGS PRIME system.
    """
    with open(fname, 'r') as f:
        x = f.readlines()
        nrows = int(x[6])
    df = pd.read_csv(fname, header=11, delimiter='\t', nrows=nrows)
    df = df.reset_index()
    df = df.rename(columns={'level_1':'a',
                            'level_3':'b',
                            'level_5':'m',
                            'level_7':'n',
                            'level_9':'resist'})
    array = df[['a','b','m','n']].values
    if espacing is None:
        espacing = np.unique(np.sort(array.flatten()))[1]
    array = np.round(array/espacing+1).astype(int)
    imax = int(np.max(array))
#                if np.sum(array == 0) > 0:
#                    print('add 1 as there is electrodes at zeros')
#                    imax = imax+1
    elec = np.zeros((imax,3))
    elec[:,0] = np.arange(0,imax)*espacing
    df['ip'] = np.nan
    df[['a','b','m','n']] = array
    
    return elec, df

# test code
#elec, df = primeParser('api/test/primeFile.dat')
    
def primeParserTab(fname, espacing = 1):
    """
    Parses data from a prime file with the .tab extension - jamyd91
    """
    #error check
    fh = open(fname,'r')
    line1 = fh.readline()
    fh.close()
    #print(line1)
    if line1.find("Prime") == -1:
        raise ImportError("Unrecognised prime file type")
        
    temp = pd.read_csv(fname,header=26,delimiter='\t')
    #Note R2 expects the electrode format in the form:
    #meas.no | P+ | P- | C+ | C- | transfer resistance
    
    a = temp["pt_p1_no:"]
    b = temp["pt_p2_no:"]
    m = temp["pt_c1_no:"]
    n = temp["pt_c2_no:"]
    num_meas = len(a)
    data_dict = {'a':a,'b':b,'m':m,'n':n}
    data_dict['resist'] = temp["pt_calc_res:"]
    data_dict['vp'] = temp["pt_meas_applied_voltage:"]
    data_dict['dev'] = temp["pt_calc_res_error:"]
    data_dict["Rho"] = [float("nan")]*num_meas
    data_dict["ip"] = [0]*num_meas
    df = pd.DataFrame(data=data_dict) # make a data frame from dictionary
    df = df[['a','b','m','n','Rho','dev','ip','resist']] # reorder columns to be consistent with the syscal parser
    
    #compute default electrode positions 
    imax = np.max(np.array((a,b,m,n)))
    elec = np.zeros((imax,3))
    elec[:,0] = np.arange(0,imax)*espacing
    return elec, df
##Test code
#electrodes, data = primeParserTab("../Data/PRIME/3001_CAN_2017-11-16_005119.tab")

#%% parse input for res2inv (.dat file) - Jimmy B. 
#jamyd91@bgs.ac.uk
def res2invInputParser(file_path):
    """
    Returns info on the electrode geometry and transfer resistances held in the res2dinv input file. 
    It looks for the general array format in the .dat file. 
    
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
    
    ## TODO : add capacity to read in borehole surveys 
    """
    fh = open(file_path,'r')#open file handle for reading
    dump = fh.readlines()#cache file contents into a list
    fh.close()#close file handle, free up resources
    
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
            print('errors estimates found in file and will be converted into transfer resistance error estimates')
        
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
        ez_pos=[0]*num_elec # actaully we can't have a z coordinate for 2d data so these will remain as zero
        for i in range(num_elec):
            ex_pos[i] = float(dump[topo_flag_idx+2+i].strip().split(',')[0])
            ey_pos[i] = float(dump[topo_flag_idx+2+i].strip().split(',')[1])
        #print(ex_pos,ey_pos)
        elec = np.column_stack((ex_pos,ey_pos,ez_pos))
              
    #since we dont always have all the electrode indexes we need to determine this
    #idea here is to extract all the x locations from the general array format
    total_x=np.array(())
    #print('finding general array electrode coordinates')
    for k in range(num_meas):
        line = dump[k+idx_oi]
        vals = line.strip().split()
        x_dump = [float(vals[1:5][i].split(',')[0]) for i in range(4)] # extract x locations
        total_x = np.append(total_x,x_dump) # this caches all the x locations 
    #extract the unique x values using numpy
    ex_pos = np.unique(total_x)
    #now we have indexed electrode coordinates in ex_pos :) 
    if not topo_flag: # then we dont have any topography and the electrode positions are simply given by thier x coordinates
        ey_pos=[0]*num_elec
        ez_pos=[0]*num_elec  
        elec = np.column_stack((ex_pos,ey_pos,ez_pos))
    
    # loke's general array format is in the form:
    #no.electrodes | C+ | C- | P+ | P- | apparent.resistivity. 
    #Note R2 expects the electrode format in the form:
    #meas.no | P+ | P- | C+ | C- | transfer resistance
    #print('computing transfer resistances and reading in electrode indexes')
    data_dict = {'a':[],'b':[],'m':[],'n':[],'Rho':[],'ip':[],'resist':[],'dev':[]}
    for k in range(num_meas):
        line = dump[k+idx_oi]
        vals = line.strip().split()
        x_dump = [float(vals[1:5][i].split(',')[0]) for i in range(4)]
        #convert the x electrode coordinates into indexes?
        e_idx = [np.where(ex_pos == x_dump[i])[0][0] for i in range(4)]
        #add the electrode indexes to the dictionary which will be turned into a dataframe
        data_dict['a'].append(e_idx[2]+1)
        data_dict['b'].append(e_idx[3]+1)
        data_dict['m'].append(e_idx[1]+1)
        data_dict['n'].append(e_idx[0]+1)
        #convert apparent resistivity back in to transfer resistance
        K = geom_fac(x_dump[0],x_dump[1],x_dump[2],x_dump[3])
        Pa = float(vals[5]) # apparent resistivity value
        Pt = Pa/K # transfer resistance
        #add apparent and transfer resistances to dictionary
        data_dict['Rho'].append(Pa)
        data_dict['resist'].append(Pt)
        data_dict['ip'].append(0)
        if err_flag:
            err_Pa = float(vals[6])
            err_Pt = err_Pa/K
            data_dict['dev'].append(abs(err_Pt))
        else:
            data_dict['dev'].append(0)
    
    df = pd.DataFrame(data=data_dict) # make a data frame from dictionary
    df = df[['a','b','m','n','Rho','dev','ip','resist']] # reorder columns to be consistent with the syscal parser
    
    return elec,df

#%% convert a dataframe output by the parsers into a simple protocal.dat file    
###DEPRECIATED###
def dataframe2dat(df,save_path='default',
                  typ='R2',
                  ref_flag=False,
                  err_flag=False):
    """
    converts a pandas dataframe into a protocol.dat file for R2.exe. 
    
    Parameters
    -----------
    df: pandas dataframe
        dataframe output by a r2gui parsers
    save_path: string
        file path to save location, if left default 'protocal.dat' is written to the working directory 
    ref_flag: boolian
        not working currently. if a reference set of transfer resistances is to be used then set to true
    err_flag: boolian
        if errors are going to be written to the r2 file. Currently uses 'dev' key in the pandas dataframe
    
    Returns
    -----------
    protocal.dat written to specificied folder/directory
    """
    num_meas = len(df)
    
    if save_path == 'default':
        save_path = 'protocol.dat'
    else:
        print(os.path.join(save_path,'protocol.dat'))
        save_path=os.path.join(save_path,'protocol.dat')
        
    #the contents of the last 3 columns will depend on survey type and regularisation mode. 
    if typ=='R2':
        column6=df['resist']
    else:
        column6=df['ip']
    
    if ref_flag:
        column7=[0]*len(column6)# placeholder -- need to identify what jkl calls his reference values 
        print("ref flag argument currently not supported and has been ignored")
    else:
        column7=[0]*len(column6)
    
    if err_flag:
        column8=df['dev']
    else:
        column8=[0]*len(column6)
     
    fh = open(save_path,'w')
    
    fh.write("%i\n"%num_meas)
    for i in range(num_meas):
        fh.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                 i+1, #measurement number 
                 df['a'][i],
                 df['b'][i],
                 df['m'][i],
                 df['n'][i],
                 column6[i],
                 column7[i],
                 column8[i]))
        
    fh.close()