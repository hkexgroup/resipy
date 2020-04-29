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

#%% function to compute geometric factor - Jamyd91
def geom_fac(C1,C2,P1,P2):
    """Compute geometric factor, as many formats give resistivities in terms of 
    apparent resistivity. R2 reads in the raw transfer resistances. 
    Parameters
    ----------
    C1: float, np array
        x position of postive current electrode
    C2: float, np array
        x position of negative current electrode
    P1: float, np array
        x position of postive potential electrode
    P2: float, np array
        x position of negative potential electrode
        
    Returns
    -----------
    k: float, np array
        geometric factor to convert transfer resistance into apparent resistivity 
    """
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
        if spacing is None:
            spacing = espacing
        array = np.round(array/spacing+1).astype(int)
        df[['a','b','m','n']] = array
        df['resist'] = df['vp']/df['i']
        imax = int(np.max(array))
        elec = np.zeros((imax,3))
        elec[:,0] = np.arange(0,imax)*spacing
                
        return elec, df
    
#test code
#elec, df = syscalParser('test/syscalFile.csv')
#elec, df = syscalParser('/media/jkl/data/phd/tmp/projects/ahdb/survey2018-08-14/data/ert/18081401.csv', 0.5)
#print(df[['a','b','m','n']])
#print(elec)

#%% protocol.dat forward modelling parser

def protocolParser(fname):
#    x = np.genfromtxt(fname, skip_header=1) # original code doesnt work for badly behaved protocol files
#    df = pd.DataFrame(x, columns=['index','a','b','m','n','resist','appResist'])
    
    with open(fname,'r') as fh:
        num_meas = int(fh.readline().strip()) # read in first line - number of measurements 
        dump = fh.readlines()
    protocol = {'index':[0]*num_meas,# protocol dictionary 
                'a':[0]*num_meas,
                'b':[0]*num_meas,
                'm':[0]*num_meas,
                'n':[0]*num_meas,
                'resist':[0]*num_meas,
                'magErr':[0]*num_meas}
#                'appResist':[0]*num_meas} 
    #determine if apparent resistivity column is present 
#    app_resis_flag = False
#    if len(dump[0])==7:
#        app_resis_flag=True
    for i, line in enumerate(dump):
        data = line.split()
        protocol['index'][i] = int(data[0])
        protocol['a'][i] = int(data[1])
        protocol['b'][i] = int(data[2])
        protocol['m'][i] = int(data[3])
        protocol['n'][i] = int(data[4])
        protocol['resist'][i] = float(data[5])
#        if app_resis_flag:
#            protocol['appResist'][i]= float(data[6])
        if len(data) == 7:
            protocol['magErr'][i] = float(data[6])
    
    if np.mean(protocol['magErr']) !=0: 
        df = pd.DataFrame(protocol) 
    else: 
        df = pd.DataFrame(protocol).drop(['magErr'], axis = 1)
    df['ip'] = np.nan
    xElec = np.arange(np.max(df[['a','b','m','n']].values))
    elec = np.zeros((len(xElec),3))
    elec[:,0] = xElec
    return elec, df

# test code
#protocolParser('api/test/protocol.dat')


def protocolParserLME(fname):
# read LME predicted errors passed back from R, only add LME predictions to df
# Michael Tso @ 190310

    with open(fname,'r') as fh:
        num_meas = int(fh.readline().strip()) # read in first line - number of measurements 
        dump = fh.readlines()
    protocol = {'index':[0]*num_meas,# protocol dictionary 
                'a':[0]*num_meas,
                'b':[0]*num_meas,
                'm':[0]*num_meas,
                'n':[0]*num_meas,
                'resist':[0]*num_meas,
                'lmeErr':[0]*num_meas}

    for i, line in enumerate(dump):
        data = line.split()
        protocol['index'][i] = int(data[0])
        protocol['a'][i] = int(data[1])
        protocol['b'][i] = int(data[2])
        protocol['m'][i] = int(data[3])
        protocol['n'][i] = int(data[4])
        protocol['resist'][i] = float(data[5])
        if len(data) == 7:
            protocol['lmeErr'][i] = float(data[6])
    lmeError = protocol['lmeErr']
    
    return lmeError

# test code
#protocolParserLME('api/test/protocol-lmeOut.dat')



#def protocolParser2(fname): # with pandas = twice slower but same time if we use np.genfromtxt()
#    colnames = np.array(['index','a','b','m','n','resist','appResist'])
##    df = pd.read_csv(fname, header=None, sep='\s+', skiprows=1)
##    colindex = 1+np.arange(df.shape[0])
##    df = df.rename(columns=dict(zip(colindex, colnames[:df.shape[0]])))
#    
#    x = np.genfromtxt(fname, skip_header=1)
#    df = pd.DataFrame(x, columns=colnames[:x.shape[1]])
#    df['ip'] = np.nan
#    xElec = np.arange(np.max(df[['a','b','m','n']].values))
#    elec = np.zeros((len(xElec),3))
#    elec[:,0] = xElec
#    return df, elec


def protocolParserIP(fname): # just for protocol forward output with cR2
#    x = np.genfromtxt(fname, skip_header=1)
#    df = pd.DataFrame(x, columns=['index','a','b','m','n','resist','ip','appResist'])
    with open(fname,'r') as fh:
        num_meas = int(fh.readline().strip()) # read in first line - number of measurements 
        dump = fh.readlines()
    protocol = {'index':[0]*num_meas,# protocol dictionary 
                'a':[0]*num_meas,
                'b':[0]*num_meas,
                'm':[0]*num_meas,
                'n':[0]*num_meas,
                'resist':[0]*num_meas,
                'ip':[0]*num_meas,
                'magErr':[0]*num_meas,
                'phiErr':[0]*num_meas}
#                'appResist':[0]*num_meas} 
    #determine if apparent resistivity column is present 
#    app_resis_flag = False
#    if len(dump[0])==8:
#        app_resis_flag=True
    for i, line in enumerate(dump):
        data = line.split()
        protocol['index'][i] = int(data[0])
        protocol['a'][i] = int(data[1])
        protocol['b'][i] = int(data[2])
        protocol['m'][i] = int(data[3])
        protocol['n'][i] = int(data[4])
        protocol['resist'][i] = float(data[5])
        protocol['ip'][i] = float(data[6])
#        if app_resis_flag:
#            protocol['appResist'][i]= float(data[7])
        if len(data) == 9:
            protocol['magErr'][i] = float(data[7])
            protocol['phiErr'][i] = float(data[8])
    
    if np.mean(protocol['phiErr']) and np.mean(protocol['magErr']) !=0: 
        df = pd.DataFrame(protocol) 
    else: 
        df = pd.DataFrame(protocol).drop(['magErr','phiErr'], axis = 1)
    xElec = np.arange(np.max(df[['a','b','m','n']].values))
    elec = np.zeros((len(xElec),3))
    elec[:,0] = xElec
    return elec, df    


#%% 3D protocol parser
def protocol3DParser(fname): # works for 2D and 3D (no IP)
    colnames3d = np.array(['index','sa', 'a','sb','b','sm', 'm','sn', 'n','resist'])
    colnames = np.array(['index','a','b','m','n','resist','magErr'])
    x = np.genfromtxt(fname, skip_header=1)
    if x.shape[1] > 8: # max for non 3D cases
        colnames = colnames3d
        if x.shape[1] > len(colnames3d): # we have IP in the file (or error)
            colnames = np.r_[colnames, ['ip']]
        # remove string effect (faster option but might lead to non-consecutive elec number)
#        maxElec = np.max(x[:,[2,4,6,8]])
#        x[:,[2,4,6,8]] = x[:,[1,3,5,7]]*maxElec + x[:,[2,4,6,8]]
#        x[:,[1,3,5,7]] = 1
#        
        # slower option (but consecutive electrode numbers)
        elecs = x[:,1:9].reshape((-1,2))
        c = 0
        for line in np.unique(elecs[:,0]):
            ie = elecs[:,0] == line
            uelec = np.unique(elecs[ie,1])
            if len(uelec) != np.max(elecs[ie,1]):
                raise Exception('More electrodes per line ({}) then electrode '
                                'number ({}).'.format(np.max(elecs[ie,1]), len(uelec)))
            else:
                elecs[ie,0] = 0
                elecs[ie,1] = c + elecs[ie,1]
                c = c + len(uelec)
        x[:,1:9] = elecs.reshape((-1,8))
        
    df = pd.DataFrame(x, columns=colnames[:x.shape[1]])
    if 'ip' not in df.columns:
        df['ip'] = np.nan
    xElec = np.arange(np.max(df[['a','b','m','n']].values))
    elec = np.zeros((len(xElec),3))
    elec[:,0] = xElec
    return elec, df

# test code
#elec1, df1 = protocol3DParser('api/test/protocol3Df.dat')
#elec2, df2 = protocol3DParser('api/test/protocol3Di.dat')
#elec3, df3 = protocol3DParser('api/test/protocol.dat')


#%% forwardProtocolDC/IP parser
    
def forwardProtocolDC(fname): # need specific as there is a appRes column
    x = np.genfromtxt(fname, skip_header=1)
    df = pd.DataFrame(x, columns=['num','a','b','m','n','resist','appRes'])
    df['ip'] = np.nan
    xElec = np.arange(np.max(df[['a','b','m','n']].values))
    elec = np.zeros((len(xElec),3))
    elec[:,0] = xElec
    return elec, df
    

def forwardProtocolIP(fname): # not needed as the normal parser can do it
    x = np.genfromtxt(fname, skip_header=1)
    df = pd.DataFrame(x, columns=['num','a','b','m','n','resist','ip','appRes'])
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

def res2InvGeneralArray(file_path):
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
    num_meas = int(dump[6])
    start = 9
    #general array format given in terms of xz coordinates 
    xC1 = np.array([0]*num_meas,dtype=float)
    xC2 = np.array([0]*num_meas,dtype=float)
    xP1 = np.array([0]*num_meas,dtype=float)
    xP2 = np.array([0]*num_meas,dtype=float)
    zC1 = np.array([0]*num_meas,dtype=float)
    zC2 = np.array([0]*num_meas,dtype=float)
    zP1 = np.array([0]*num_meas,dtype=float)
    zP2 = np.array([0]*num_meas,dtype=float)
    pa = np.array([0]*num_meas,dtype=float)
    Tr = np.array([0]*num_meas,dtype=float)
    count = 0
    for i in range(start,num_meas+start):
        line = dump[i].split()
        xC1[count] = float(line[1])
        zC1[count] = float(line[2])
        xC2[count] = float(line[3])
        zC2[count] = float(line[4])
        xP1[count] = float(line[5])
        zP1[count] = float(line[6])
        xP2[count] = float(line[7])
        zP2[count] = float(line[8])
        pa[count] = float(line[9])
        try:
            k = geom_fac(xC1[count],xC2[count],xP1[count],xP2[count])
        except:
            k=float('nan')
        Tr[count] = pa[count]/k
        count += 1
        
    #convert xz coordinates into electrode numbers 
    uni_x, inds = np.unique(np.column_stack((xC1,xC2,xP1,xP2)).flatten(),return_index=True)
    total_z = np.column_stack((zC1,zC2,zP1,zP2)).flatten()
    uni_z = total_z[inds]
    uni_list = list(uni_x)
    num_elec = len(uni_x)
    
    C1 = np.array([0]*num_meas)
    C2 = np.array([0]*num_meas)
    P1 = np.array([0]*num_meas)
    P2 = np.array([0]*num_meas)
    for i in range(num_meas):
        C1[i] = uni_list.index(xC1[i])+1
        C2[i] = uni_list.index(xC2[i])+1
        P1[i] = uni_list.index(xP1[i])+1
        P2[i] = uni_list.index(xP2[i])+1
    
    elec = np.zeros((num_elec,3))
    elec[:,0] = uni_x
    elec[:,1] = uni_z
    
    #return output in dataframe format
    data_dict = {'a':[],'b':[],'m':[],'n':[],'Rho':[],'ip':[],'resist':[],'dev':[]}
    data_dict['a']=P1
    data_dict['b']=P2
    data_dict['n']=C2
    data_dict['m']=C1
    data_dict['resist']=Tr
    data_dict['Rho']=pa
    data_dict['dev']=[0]*num_meas
    data_dict['ip']=[0]*num_meas
    df = pd.DataFrame(data=data_dict) # make a data frame from dictionary
    df = df[['a','b','m','n','Rho','dev','ip','resist']] # reorder columns to be consistent with the syscal parser
    return elec,df
    
#return dipole dipole survey 
def res2InvDipoleDipole(fname):
    """Read in 2D ABEM file for a dipole dipole array ONLY. 
    """
    fh = open(fname,'r')
    dump = fh.readlines()
    fh.close()
    
    num_meas = int(dump[3])#number of measurements 
    a_factor = [0]*num_meas#a spacing value
    n_factor = [0]*num_meas#n spacing value
    xpos = [0]*num_meas#leftmost electrode x position 
    pa = [0]*num_meas# apparent resisitivity 
    
    for i in range(num_meas):
        line = dump[6+i].split()
        xpos[i] = float(line[0])
        a_factor[i] = float(line[1])
        n_factor[i] = float(line[2])
        pa[i] = float(line[3])
     
    end_meas = 6+num_meas
    num_elec = int(dump[end_meas+1])
    
    #get elec coordinates 
    elec_x =  [0]*num_elec
    elec_z =  [0]*num_elec   
    elec_id =  [0]*num_elec
    
    for i in range(num_elec):
        line = dump[end_meas+2+i].split()
        elec_x[i] = float(line[0])
        elec_z[i] = float(line[1])
        elec_id[i] = i+1
        
    #elec = np.array([elec_x,[0]*num_elec,elec_z]).T
    elec = np.array([elec_x,elec_z,[0]*num_elec]).T
    # convert to quadrapoles 
        
    # ABEM general array format is in the form:
    #no.electrodes | C- | C+ | P+ | P- | apparent.resistivity. 
    #Note R2 expects the electrode format in the form:
    #meas.no | P+ | P- | C+ | C- | transfer resistance
    
    #a,b,n,m
    a=[0]*num_meas
    b=[0]*num_meas
    n=[0]*num_meas
    m=[0]*num_meas
    Tr =[0]*num_meas # transfer resistance
    
    for i in range(num_meas):
        c2_loc = xpos[i] # C- 
        c1_loc = xpos[i]+a_factor[i] # C+
        p1_loc = c1_loc + (a_factor[i]*n_factor[i]) # p+
        p2_loc = p1_loc + a_factor[i] #P- 
        a[i] = elec_x.index(p1_loc)+1
        b[i] = elec_x.index(p2_loc)+1
        m[i] = elec_x.index(c2_loc)+1
        n[i] = elec_x.index(c1_loc)+1
        K = geom_fac(c1_loc,c2_loc,p1_loc,p2_loc)
        Tr[i] = pa[i]/K
    
    #put data into correct format
    data_dict = {'a':[],'b':[],'m':[],'n':[],'Rho':[],'ip':[],'resist':[],'dev':[]}
    data_dict['a']=a
    data_dict['b']=b
    data_dict['n']=n
    data_dict['m']=m
    data_dict['resist']=Tr
    data_dict['Rho']=pa
    data_dict['dev']=[0]*num_meas
    data_dict['ip']=[0]*num_meas
    df = pd.DataFrame(data=data_dict) # make a data frame from dictionary
    df = df[['a','b','m','n','Rho','dev','ip','resist']] # reorder columns to be consistent with the syscal parser
    
    return elec,df

def resInvParser(filename):
    """Returns info on the electrode geometry and transfer resistances held in the res2dinv input file. 
    
    Parameters
    -----------
    file_path : string 
         string mapping to the res2inv input file 
    """
    fh = open(filename,'r')
    _ = fh.readline()#title line
    _ = fh.readline()
    flag = int(fh.readline())
    fh.close()
    if flag == 3: 
        elec, df = res2InvDipoleDipole(filename)
    elif flag == 11:
        elec, df = res2InvGeneralArray(filename)
    elif flag != 3 or flag != 11:
        elec, df = res2invInputParser(filename)
    else:
        raise ImportError("unsupported type of res inv file (for the moment)")
        
    return elec,df

#%% parse 3D sting data

def stingParser(fname):
    """Read in .stg file from sting
    """
    fh = open(fname,'r')
    dump = fh.readlines()
    fh.close()    
    
    num_meas = int(dump[1].split()[-1]) # number of measurements 
    Tr = [0]*num_meas # transfer resistance - pre allocate arrays to populate after read in 
    pa = [0]*num_meas # apparent resistivity 
    meas_id = [0]*num_meas # measurement id number
    c1_loc = [0]*num_meas
    c2_loc = [0]*num_meas
    p1_loc = [0]*num_meas
    p2_loc = [0]*num_meas
    
    for i in range(num_meas):
        line = dump[i+3].split(',')
        Tr[i] = float(line[4])
        c1_loc[i] = float(line[12])#C+
        c2_loc[i] = float(line[9])#C-
        p1_loc[i] = float(line[15])#P+
        p2_loc[i] = float(line[18])#P-
        meas_id[i] = int(line[0])
        pa[i] = float(line[7])
    # sting array format is in the form:
    #no.electrodes | C- | C+ | P+ | P- | apparent.resistivity. 
    #Note R2 expects the electrode format in the form:
    #meas.no | P+ | P- | C+ | C- | transfer resistance
        
    loc_array=np.array([c1_loc,c2_loc,p1_loc,p2_loc]).T
    elec_x=list(np.unique(loc_array.flatten()))
    
    elec =np.array([elec_x,np.zeros_like(elec_x),np.zeros_like(elec_x)]).T # electrode array
    #TODO: return true positions? 
    
    a = [0]*num_meas
    b = [0]*num_meas
    n = [0]*num_meas
    m = [0]*num_meas
    
    for i in range(num_meas):
        a[i] = elec_x.index(p1_loc[i])+1
        b[i] = elec_x.index(p2_loc[i])+1
        m[i] = elec_x.index(c2_loc[i])+1
        n[i] = elec_x.index(c1_loc[i])+1
    
    #put data into correct format
    data_dict = {'a':[],'b':[],'m':[],'n':[],'Rho':[],'ip':[],'resist':[],'dev':[]}
    data_dict['a']=a
    data_dict['b']=b
    data_dict['n']=n
    data_dict['m']=m
    data_dict['resist']=Tr
    data_dict['Rho']=pa
    data_dict['dev']=[0]*num_meas
    data_dict['ip']=[0]*num_meas
    df = pd.DataFrame(data=data_dict) # make a data frame from dictionary
    df = df[['a','b','m','n','Rho','dev','ip','resist']] # reorder columns to be consistent with the syscal parser
    
    return elec,df

#fname = '070708L5_trial1.stg'
#stingParser(fname)
