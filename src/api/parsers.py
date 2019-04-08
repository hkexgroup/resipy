#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  1 11:23:23 2018
Handles importing external data files into the R2 api
@authors: Guillaume, Jimmy, Sina, Pedro Concha 

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
    Rc1p1 = np.abs(C1 - P1)
    Rc2p1 = np.abs(C2 - P1)
    Rc1p2 = np.abs(C1 - P2)
    Rc2p2 = np.abs(C2 - P2)
    
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
    
    ## TODO : add capacity to read Pole-Pole, Pole-Dipole and borehole surveys 
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
    
    meas_type_flag = 'appRes' #default  
    if array_type in [1,3,4,5,7]: 
        #1:Wenner alfa, 3:DipDip,4:Wenner beta, 5:Wenner gama,
        #7:Schlumberger, 11:General array, 15:Gradient 
        idx_oi = 3
        line = dump[idx_oi]
        
    elif array_type in (11, 15):
        idx_oi = 6
        line = dump[idx_oi]
    
        meas_type = [int(dump[i+1]) for i in range(len(dump)) if 'Type' in dump[i]] # Looking for app.resistivity/resistiance flag
           
        if meas_type != []:
            if meas_type[0] == 1:
                meas_type_flag = 'resistance'

        elif meas_type == []:
            if dump[4] == 1 or dump[5] == 1:
                meas_type_flag = 'resistance'
    
    num_meas = int(line)
    idx_oi += 1
    line = dump[idx_oi]
    x_location = int(line)
    idx_oi += 1

    ip_flag = False
    for i in range(len(dump)): 
        if 'Chargeability' in dump[i]:
            ip_flag = True
            factor = 1
        elif 'Phase' in dump[i]:
            ip_flag = True
            factor = -1
            
    if ip_flag == True:
        idx_oi += 4
    else:
        idx_oi += 1
    
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
            c2 = (float(vals[3]))
            p1 = (float(vals[5]))
            p2 = (float(vals[7]))
            
        x_dump.append(c1)
        x_dump.append(p1)
        x_dump.append(p2)
        x_dump.append(c2)
        total_x = np.append(total_x,x_dump)
        total_x = np.around(total_x, decimals = 1)
        #convert the x electrode coordinates into indexes?
        ex_pos = np.unique(total_x)
        ex_pos = np.around(ex_pos, decimals = 1)
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

        elif array_type in (11,15):
            c1 = (float(vals[1]))
            c2 = (float(vals[3]))
            p1 = (float(vals[5]))
            p2 = (float(vals[7]))
            Pa.append(float(vals[9]))
        if any([c1 == 0,c2 == 0,p1 == 0,p2 == 0]):
            start_0_flag = True
        c1 = round(c1, 0)
        p1 = round(p1, 0)
        p2 = round(p2, 0)
        c2 = round(c2, 0)
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
            K = geom_fac(c1, c2, p1, p2)
        elif array_type == 3:
            K = np.pi * n*(n + 1)*(n + 2)*a
        elif array_type == 4:
            K = 6 * np.pi * a
        elif array_type == 5:
            K = 3 * np.pi * a

        #add apparent and transfer resistances to dictionary
        if meas_type_flag == 'resistance':
            R = Pa[k] # transfer resistance
        else:
            R = Pa[k]/K
            
        data_dict['Rho'].append(Pa[k])
        data_dict['resist'].append(R)
        if ip_flag == True:
            data_dict['ip'].append(factor*float(vals[-1])) # add ip (chargeability) assuming it's last column
        else:
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


def resInvParser(filename): # keeping it for now, in case of Res3DInv files
    """Returns info on the electrode geometry and transfer resistances held in the res2dinv input file. 
    
    Parameters
    -----------
    filename : string 
         string mapping to the res2inv input file 
    """
    try:
        elec, df = res2invInputParser(filename)
    except:
        raise ImportError('Unsupported ResInv file')
        
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
