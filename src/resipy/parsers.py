#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  1 11:23:23 2018
Parse device-specific format and extract data as pandas.Dataframe 
and electrode positions as numpy.array
@authors: Guillaume, Jimmy, Sina, Pedro Concha
"""

import numpy as np
import pandas as pd
import os 
import re

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

#%% Functions for ericParser
def ndmesh(*xi,**kwargs):
    if len(xi) < 2:
        msg = 'meshgrid() takes 2 or more arguments (%d given)' % int(len(xi) > 0)
        raise ValueError(msg)

    args = np.atleast_1d(*xi)
    ndim = len(args)
    copy_ = kwargs.get('copy', True)

    s0 = (1,) * ndim
    output = [x.reshape(s0[:i] + (-1,) + s0[i + 1::]) for i, x in enumerate(args)]

    shape = [x.size for x in output]

    # Return the full N-D matrix (not only the 1-D vector)
    if copy_:
        mult_fact = np.ones(shape, dtype=int)
        return [x * mult_fact for x in output]
    else:
        return np.broadcast_arrays(*output)
    
    
def geom_factor_3D(df, elec, array_type):
    """Compute geometric factor
    """
    array = df[['a','b','m','n']].values.astype(int)
#    array = df[['a','b','m','n']].values
    aposx = elec[:,0][array[:,0]-1]
    aposy = elec[:,1][array[:,0]-1]
    aposz = elec[:,2][array[:,0]-1]
    if array_type != [4]:    
        bposx = elec[:,0][array[:,1]-1]
        bposy = elec[:,1][array[:,1]-1]
        bposz = elec[:,2][array[:,1]-1]
        
    mposx = elec[:,0][array[:,2]-1]
    mposy = elec[:,1][array[:,2]-1]
    mposz = elec[:,2][array[:,2]-1]
    if array_type != [4]:     
        nposx = elec[:,0][array[:,3]-1]
        nposy = elec[:,1][array[:,3]-1]
        nposz = elec[:,2][array[:,3]-1]
        
    am = np.sqrt((aposx-mposx)**2 + (aposy-mposy)**2 + (aposz-mposz)**2)
    if array_type != [4]: 
        bm = np.sqrt((bposx-mposx)**2 + (bposy-mposy)**2 + (bposz-mposz)**2)
        an = np.sqrt((aposx-nposx)**2 + (aposy-nposy)**2 + (aposz-nposz)**2)
        bn = np.sqrt((bposx-nposx)**2 + (bposy-nposy)**2 + (bposz-nposz)**2)
    if array_type == [4]:
        k = 2*np.pi/(1/am)
    else:
        k = 2*np.pi/((1/am)-(1/bm)-(1/an)+(1/bn)) # geometric factor
                
    return k

#%% usual syscal parser
def syscalParser(fname):#, spacing=None):
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
    
        df['resist'] = df['vp']/df['i']

        # if spacing is None:    
        # for unregularly spaced array
        array = df[['a','b','m','n']].values
        # arrayMin = np.min(np.unique(np.sort(array.flatten())))
        # if arrayMin != 0: # all surveys must start from x = 0
        #     array -= arrayMin
        # val = np.sort(np.unique(array.flatten())) # unique electrodes positions
        # elecLabel = 1 + np.arange(len(val))
        # newval = elecLabel[np.searchsorted(val, array)] # magic ! https://stackoverflow.com/questions/47171356/replace-values-in-numpy-array-based-on-dictionary-and-avoid-overlap-between-new
        # df.loc[:,['a','b','m','n']] = newval
        # elec = np.c_[val, np.zeros((len(val),2))]
        remoteFlags = np.array([-9999999, -999999, -99999,-9999,-999,
                                9999999, 999999, 99999, 9999, 999])
        iremote = np.in1d(array.flatten(), remoteFlags)
        remoteFreeArray = array.flatten()[~iremote]
        arraySorted = np.sort(np.unique(remoteFreeArray))
        elecSpacing = arraySorted[1] - arraySorted[0]
        # print('elecSpacing=', elecSpacing)
        n = np.max(remoteFreeArray/elecSpacing).astype(int)
        a = 0
        if arraySorted[0] == 0:
            n = n + 1
            a = 1
        if np.sum(iremote) > 0:
            n = n + 1
        for i in range(4):
            ie = ~np.in1d(array[:,i], remoteFlags)
            array[ie, i] = array[ie, i]/elecSpacing + a
            array[~ie, i] = n
        df.loc[:,['a','b','m','n']] = array.astype(int)
        elec = np.zeros((n, 3))
        elec[:,0] = np.arange(elec.shape[0])*elecSpacing # assumed regular spacing                
        if np.sum(iremote) > 0:
            elec[-1, 0] = -99999
        # else:
        #     # for regularly spaced array
        #     array = df[['a','b','m','n']].values
        #     arrayMin = np.min(np.unique(np.sort(array.flatten())))
        #     if arrayMin != 0:
        #         array -= arrayMin
        #     espacing = np.unique(np.sort(array.flatten()))[1] - np.unique(np.sort(array.flatten()))[0]
        #     if spacing is None:
        #         spacing = espacing
        #     array = np.round(array/spacing+1).astype(int)
        #     df[['a','b','m','n']] = array
        #     imax = int(np.max(array))
        #     elec = np.zeros((imax,3))
        #     elec[:,0] = np.arange(0,imax)*spacing
                
        return elec, df
    
#test code
# elec, df = syscalParser('examples/dc-2d/syscal.csv')
# elec, df = syscalParser('examples/dc-2d-pole-dipole/syscal.csv')

#%%
def protocolParserLME(fname): # pragma: no cover
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



#%% protocol parser for 2D/3D and DC/IP
def protocolParser(fname, ip=False, fwd=False):
    """
    <type>     <ncol>
    DC 2D         6
    DC 2D + err   7
    DC 2D + fwd   7
    IP 2D         7
    IP 2D + err   9
    IP 2D + fwd   8
    DC 3D         10
    DC 3D + err   11
    DC 3D + fwd   11
    IP 3D         11
    IP 3D + err   13
    IP 3D + fwd   12
    """
    x = np.genfromtxt(fname, skip_header=1) # we don't know if it's tab or white-space
    if fwd:
        x = x[:,:-1] # discard last column as it is appRes
    if ip:
        colnames3d = np.array(['index','sa','a','sb','b','sm', 'm','sn','n','resist','ip','magErr','phiErr'])
        colnames2d = np.array(['index','a','b','m','n','resist','ip','magErr','phiErr'])
    else:
        colnames3d = np.array(['index','sa','a','sb','b','sm', 'm','sn','n','resist','magErr'])
        colnames2d = np.array(['index','a','b','m','n','resist','magErr'])
    ncol = x.shape[1]
    if ncol <= len(colnames2d): # it's a 2D survey
        colnames = colnames2d[:ncol]
    else: # it's a 3D survey
        colnames = colnames3d[:ncol]
        # putting all electrodes on the same string            
        # lineNum = x[:,1:9:2] # line numbers 
        # elecNum = x[:,2:9:2] # electrode numbers 
        # lineNumF = lineNum.flatten() # flattened arrays 
        # elecNumF = elecNum.flatten()
        # c = 0
        # for line in np.unique(lineNumF):
        #     ie = lineNumF == line # electrode indexes 
        #     elecNumF[ie] += c # add maximum electrode index found so far 
        #     c = np.max(elecNumF[ie])
        # measNum = x.shape[0] # number of measurements 
        # x[:,1:9:2] = np.ones((measNum,4)) # make line numbers all 1
        # x[:,2:9:2] = elecNumF.reshape(elecNum.shape)
        
    df = pd.DataFrame(x, columns=colnames)
    df = df.astype({'a':int, 'b':int, 'm':int, 'n':int})
    if 'sa' in df.columns:
        df = df.astype({'sa':int, 'sb':int, 'sm':int, 'sn':int})
        elec = np.vstack([df[['sa','a']].values, df[['sb','b']].values,
                          df[['sm','m']].values, df[['sn','n']].values])
        uelec = np.unique(elec, axis=0)
        dfelec = pd.DataFrame(uelec, columns=['string', 'elec'])
        dfelec = dfelec.sort_values(by=['string','elec']).reset_index(drop=True)
        dfelec['label'] = dfelec['string'].astype(str) + ' ' + dfelec['elec'].astype(str)
        dfelec = dfelec.drop(['string', 'elec'], axis=1)
        df['a2'] = df['sa'].astype(str) + ' ' + df['a'].astype(str)
        df['b2'] = df['sb'].astype(str) + ' ' + df['b'].astype(str)
        df['m2'] = df['sm'].astype(str) + ' ' + df['m'].astype(str)
        df['n2'] = df['sn'].astype(str) + ' ' + df['n'].astype(str)
        df = df.drop(['a','b','m','n','sa','sb','sm','sn'], axis=1)
        df = df.rename(columns={'a2':'a','b2':'b','m2':'m','n2':'n'})
    else:
        uelec = np.unique(df[['a','b','m','n']].values.flatten()).astype(int)
        dfelec = pd.DataFrame(uelec, columns=['label'])
        dfelec = dfelec.astype({'label': str})
    dfelec['x'] = np.arange(dfelec.shape[0])
    dfelec['y'] = 0
    dfelec['z'] = 0
    dfelec['buried'] = False
    dfelec['remote'] = False
    df = df.astype({'a':str, 'b':str, 'm':str, 'n':str})
    if 'ip' not in df.columns:
        df['ip'] = np.nan
    return dfelec, df

# test code
# elec, df = protocolParser('examples/dc-2d/protocol.dat')
# elec, df = protocolParser('examples/dc-3d/protocol.dat')
# elec, df = protocolParser('examples/ip-2d/protocol.dat', ip=True)
# elec, df = protocolParser('examples/ip-3d/protocol2.dat', ip=True)


# s1 = np.unique(elec['label'].values)
# s2 = np.unique(df[['a','b','m','n']].values.flatten())
# x = np.intersect1d(s1, s2)
# lookupDict = dict(zip(elec['label'], np.arange(elec.shape[0])))
# array = df[['a','b','m','n']].replace(lookupDict).values
# print(array)


#%% PRIME system parser

def primeParserTab(fname, espacing=1):
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
    
    # shift electrodes to have continuous numbering
    # array = df[['a','b','m','n']].values
    # arrayMin = np.min(np.unique(np.sort(array.flatten())))
    # if arrayMin != 0: # all surveys must start from x = 0
    #     array -= arrayMin
    # val = np.sort(np.unique(array.flatten())) # unique electrodes positions
    # elecLabel = 1 + np.arange(len(val))
    # newval = elecLabel[np.searchsorted(val, array)] # magic ! https://stackoverflow.com/questions/47171356/replace-values-in-numpy-array-based-on-dictionary-and-avoid-overlap-between-new
    # df.loc[:,['a','b','m','n']] = newval
    
    #compute default electrode positions
    array = df[['a','b','m','n']].values
    arraySorted = np.sort(np.unique(array.flatten()))
    n = np.max(arraySorted)
    if arraySorted[0] == 0:
        n = n + 1
    elec = np.zeros((n, 3))
    elec[:,0] = np.arange(n) * espacing
    
    return elec, df
##Test code
#electrodes, data = primeParserTab("../Data/PRIME/3001_CAN_2017-11-16_005119.tab")


#%% parse input for res2inv (.dat file)
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
    
    ## TODO : add capacity to read Offset Pole-Dipole and borehole surveys 
    """
    c1 = np.array(())
    c2 = np.array(())
    p1 = np.array(())
    p2 = np.array(())
    c1_y = np.array(())
    c2_y = np.array(())
    p1_y = np.array(())
    p2_y = np.array(())
    c1_z = np.array(()) # for general array integrated topography
    c2_z = np.array(())
    p1_z = np.array(())
    p2_z = np.array(())
    pa = np.array(())
    ip = np.array(())
    r1 = np.array(()) # For approx geom factor
    r2 = np.array(())

    
    fh = open(file_path,'r')#open file handle for reading
    dump = fh.readlines()#cache file contents into a list
    fh.close()#close file handle, free up resources
    
    #first find loke's General array format information in the file (which is all we need for R2/3)
    fmt_flag = False # format flag 
    err_flag = False # error estimates flag 
    topo_flag = False # topography flag
    topo_flag_GA = False # topography flag for integrated topography in general arrays
    start_0_flag = False
    sur_flag = 0 # survey type
    x_location = 0 # for data points, 0 for first electrode, 1 for mid-point, 2 for surface distance
    ip_flag = False # 0 for none IP data, 1 if present
    factor_used = False #Exact Geometric factor used is true 

    idx_oi = 0
    line = dump[idx_oi]
    sur_name = line.strip() #name of survey
    idx_oi += 1
    line = dump[idx_oi]
    a_spac = float(line) #electrode spacing
    idx_oi += 1
    line = dump[idx_oi]
    array_type = int(line)
    #1:Wenner alfa,2:Pole_pole 3:Dipole-Dipole,4:Wenner beta, 5:Wenner gama,
    #6:Pole-Dipole 7:Schlumberger, 11:General array, 15:Gradient 
    meas_type_flag = 'appRes' #default  

    if array_type in [8,12,13]:
        raise ImportError("Not supported")
        # 8: Equatorial Dipole Dipole
        #12: Cross-borehole survey(apparent resistivity values)
        #13: Cross-borehole survey(resistance values)
   
    if array_type in [1,2,3,4,5,6,7,10]: 
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
    
    if array_type in (2, 6):
        if 'Remote electrodes included' in dump[idx_oi]:
            line = dump[idx_oi+2]
            vals = line.strip().split(',')
            c2_x_pos = float(vals[0])
            c2_y_pos = float(vals[1])
            c2_z_pos = float(vals[2])
            if array_type == 2:
                line = dump[idx_oi+4]
                vals = line.strip().split(',')
                p2_x_pos = float(vals[0])
                p2_y_pos = float(vals[1])
                p2_z_pos = float(vals[2])
            if array_type == 6:
                idx_oi = 6
            if array_type == 2:
                idx_oi = 8
            if 'Exact Geometric factor used' in dump[idx_oi]:
                factor_used = True
            else:
                factor_used = False
            idx_oi += 1
            line = dump[idx_oi]
            
    if array_type == 8 or array_type == 10:
        b_dist = float(line)
        idx_oi += 1
        line = dump[idx_oi]
            
    num_meas = int(line)
    idx_oi += 1
    line = dump[idx_oi]
    x_location = int(line)
    #x_location = 0 : First electrode location
    #x_location = 1 : Mid-point location
    #x_location = 2 : Surface distance

    idx_oi += 1
    line = dump[idx_oi]
    if int(line)== 0:
        ip_flag = False
        idx_oi += 1
    else:
        ip_flag = True
        idx_oi += 4

    total_x = np.array(())
    e_idx = []
#    total_x = np.append(total_x,0)
    data_dict = {'a':[],'b':[],'m':[],'n':[],'Rho':[],'ip':[],'resist':[],'dev':[]}
    
    for k in range(num_meas):
        line = dump[idx_oi + k]
        vals = re.findall(r'[-+]?\d*\.\d+|\d+', line)
        a = float(vals[1])
        if array_type == 1:# Wenner alfa
            pa = np.append(pa, float(vals[2]))
            if ip_flag:
                ip = np.append(ip, float(vals[3]))
            if x_location == 0:
                c1 = np.append(c1, np.around((float(vals[0])), 2))
                p1 = np.append(p1, np.around((float(vals[0]) + a), 2))
                p2 = np.append(p2, np.around((float(vals[0]) + 2*a), 2))
                c2 = np.append(c2, np.around((float(vals[0]) + 3*a), 2))
            if x_location == 1:
                mid_point = (float(vals[0]))
                p1 = np.append(p1, np.around((mid_point - a/2), 2))
                c1 = np.append(c1, np.around((mid_point - 3*a/2), 2))
                p2 = np.append(p2, np.around((mid_point + a/2), 2))
                c2 = np.append(c2, np.around((mid_point + 3*a/2), 2))
        elif array_type == 2:
            pa = np.append(pa, float(vals[2]))
            r1 = np.append(r1, float(vals[1]))
            if ip_flag:
                ip = np.append(ip, float(vals[3]))
            if x_location == 0:
                c1 = np.append(c1, float(vals[0]))
                p1 = np.append(p1, float(vals[0] + a))
                if factor_used:
                    c2_x_len = c2_x_pos - float(vals[0])
                    c2_y_len = c2_y_pos - 0
                    p2_x_len = p2_x_pos - (float(vals[0]) + a)
                    p2_y_len = p2_y_pos - 0
            if x_location == 1:
                mid_point = float(vals[0])
                c1 = np.append(c1, mid_point - a/2)
                p1 = np.append(p1, mid_point + a/2)
                if factor_used:
                    c2_x_len = c2_x_pos - (mid_point - a/2)
                    c2_y_len = c2_y_pos - 0
                    p2_x_len = p2_x_pos - (mid_point + a/2)
                    p2_y_len = p2_y_pos - 0
            if factor_used:
                c2_dist = np.sqrt(c2_x_len**2 + c2_y_len**2)
                c2 = np.append(c2, c2_dist)
                p2_dist = np.sqrt(p2_x_len**2 + p2_y_len**2)
                p2 = np.append(p2, p2_dist)
            else:
                c2 = np.append(c2,-999999)
                p2 = np.append(p2, 999999)
        elif array_type == 3:#Dipole-Dipole
            n = (float(vals[2]))
            pa = np.append(pa, float(vals[3]))
            if ip_flag:
                ip = np.append(ip, float(vals[4]))
            if x_location == 0:
                c2 = np.append(c2, np.around((float(vals[0])), 2))
                c1 = np.append(c1, np.around((float(vals[0]) + a), 2))
                p1 = np.append(p1, np.around((float(vals[0]) + a*(1 + n)), 2))
                p2 = np.append(p2, np.around((float(vals[0]) + a*(2 + n)), 2))
            if x_location == 1:
                mid_point = (float(vals[0]))
                c1 = np.append(c1, mid_point - n*a/2)
                c2 = np.append(c2, mid_point - a*((n/2) + 1))
                p1 = np.append(p1, mid_point + n*a/2)
                p2 = np.append(p2, mid_point + a*((n/2) + 1))
        elif array_type == 4:# Wenner beta
            pa = np.append(pa, float(vals[2]))
            if ip_flag:
                ip = np.append(ip, float(vals[3]))
            if x_location == 0:
                c2 = np.append(c2, np.around((float(vals[0])), 2))
                c1 = np.append(c1, np.around((float(vals[0]) + a), 2)) 
                p1 = np.append(p1, np.around((float(vals[0]) + 2*a), 2))
                p2 = np.append(p2, np.around((float(vals[0]) + 3*a), 2))
            if x_location == 1:
                mid_point = (float(vals[0]))
                c1 = np.append(c1, np.around((mid_point - a/2), 2))
                c2 = np.append(c2, np.around((mid_point - 3*a/2), 2))
                p1 = np.append(p1, np.around((mid_point + a/2), 2))
                p2 = np.append(p2, np.around((mid_point + 3*a/2), 2))
        elif array_type == 5:# Wenner gamma
            pa = np.append(pa, float(vals[2]))
            if ip_flag:
                ip = np.append(ip, float(vals[3]))
            if x_location == 0:
                c1 = np.append(c1, np.around((float(vals[0])), 2))
                p1 = np.append(p1, np.around((float(vals[0]) + a), 2))
                c2 = np.append(c2, np.around((float(vals[0]) + 2*a), 2))
                p2 = np.append(p2, np.around((float(vals[0]) + 3*a), 2))
            if x_location == 1:
                mid_point = (float(vals[0]))
                p1 = np.append(p1, np.around((mid_point - a/2), 2))
                c1 = np.append(c1, np.around((mid_point - 3*a/2), 2))
                c2 = np.append(c2, np.around((mid_point + a/2), 2))
                p2 = np.append(p2, np.around((mid_point + 3*a/2), 2))
        elif array_type == 6:# Pole-Dipole
            pa = np.append(pa, float(vals[3]))
            if ip_flag:
                ip = np.append(ip, float(vals[4]))
            n = (float(vals[2]))
            if x_location == 0:
                if n > 0:
                    c1 = np.append(c1, np.around((float(vals[0])), 2))
                    p1 = np.append(p1, np.around((float(vals[0]) + n * a), 2))
                    p2 = np.append(p2, np.around((float(vals[0]) + a*(1 + n)), 2))
                else:
                    p2 = np.append(p2, np.around((float(vals[0])), 2))
                    p1 = np.append(p1, np.around((float(vals[0]) + a), 2))
                    c1 = np.append(c1, np.around((float(vals[0]) + a*(1 + abs(n))), 2))
                if factor_used:
                    c2_x_len = c2_x_pos - float(vals[0])
                    c2_y_len = c2_x_pos - 0
            if x_location == 1:
                mid_point = (float(vals[0]))
                if n > 0:
                    c1 = np.append(c1, np.around((mid_point - n*a/2), 2))
                    p1 = np.append(p1, np.around((mid_point + n*a/2), 2))
                    p2 = np.append(p2, np.around((mid_point + a*((n/2) + 1)), 2))
                else:
                    p1 = np.append(p1, np.around((mid_point - abs(n)*a/2), 2))
                    p2 = np.append(p2, np.around((mid_point - a*(1 + (abs(n)/2))), 2))
                    c1 = np.append(c1, np.around((mid_point + abs(n)*a/2), 2))
                if factor_used:
                    c2_x_len = c2_x_pos - (np.around((mid_point + abs(n)*a/2), 2))
                    c2_y_len = c2_y_pos - 0
            if factor_used:
                c2_dist = np.sqrt(c2_x_len**2 + c2_y_len**2)
                c2 = np.append(c2, c2_dist)
            else:
                c2 = np.append(c2,-999999)
        elif array_type == 7:#Schlumberger
            pa = np.append(pa, float(vals[3]))
            if ip_flag:
                ip = np.append(ip, float(vals[4]))
            n = (float(vals[2]))
            if x_location == 0:
                c1 = np.append(c1, np.around(float(vals[0]), 2))
                p1 = np.append(p1, np.around((float(vals[0]) + (n * a)), 2))
                p2 = np.append(p2, np.around((float(vals[0]) + (a*(1 + n))), 2))
                c2 = np.append(c2, np.around((float(vals[0]) + (a*(1 + 2*n))), 2))
            if x_location == 1:
                mid_point = (float(vals[0]))
                p1 = np.append(p1, np.around((mid_point - (a/2)), 2))
                c1 = np.append(c1, np.around((mid_point - (a*(n + 1/2))), 2))
                p2 = np.append(p2, np.around((mid_point + (a/2)), 2))
                c2 = np.append(c2, np.around((mid_point + (a*(n + 1/2))), 2))
        elif array_type == 8:#Equatorial dipole-dipole >> weird!
            pa = np.append(pa, float(vals[2]))
            if ip_flag:
                ip = np.append(ip, float(vals[3]))
            if x_location == 0:
                c1 = np.append(c1, np.around((float(vals[0])), 2))
                c1_y = np.append(c1_y, 0)
                p1 = np.append(p1, np.around((float(vals[0]) + a), 2))
                p1_y = np.append(p1_y, 0)
                p2 = np.append(p2, np.around((float(vals[0]) + a), 2))
                p2_y = np.append(p2_y, b_dist)
                c2 = np.append(c2, np.around((float(vals[0])), 2))
                c2_y = np.append(c2_y, b_dist)
            if x_location == 1:
                mid_point = (float(vals[0]))
                c1 = np.append(c1, np.around((mid_point - a/2), 2))
                c1_y = np.append(c1_y, 0)
                c2 = np.append(c2, np.around((mid_point - a/2), 2))
                c2_y = np.append(c2_y, b_dist)
                p1 = np.append(p1, np.around((mid_point + a/2), 2))
                p1_y = np.append(p1_y, 0)
                p2 = np.append(p2, np.around((mid_point + a/2), 2))
                p2_y = np.append(p2_y, b_dist)
        elif array_type == 10:#Offset pole-dipole
            pa = np.append(pa, float(vals[3]))
            if ip_flag:
                ip = np.append(ip, float(vals[4]))
            n = (float(vals[2]))
            if x_location == 0:
                if n > 0:
                    c1 = np.append(c1, np.around((float(vals[0])), 2))
                    p1 = np.append(p1, np.around((float(vals[0]) + n * a), 2))
                    p2 = np.append(p2, np.around((float(vals[0]) + a*(1 + n)), 2))
                else:
                    p2 = np.append(p2, np.around((float(vals[0])), 2))
                    p1 = np.append(p1, np.around((float(vals[0]) + a), 2))
                    c1 = np.append(c1, np.around((float(vals[0]) + a*(1 + abs(n))), 2))
            if x_location == 1:
                mid_point = (float(vals[0]))
                if n > 0:
                    c1 = np.append(c1, np.around((mid_point - (a*(n + 1)/2)), 2))
                    p1 = np.append(p1, np.around((mid_point + (a*(n - 1)/2)), 2))
                    p2 = np.append(p2, np.around((mid_point + (a*(n + 1)/2)), 2))
                else:
                    p1 = np.append(p1, np.around((mid_point - (a*(abs(n) - 1)/2)), 2))
                    p2 = np.append(p2, np.around((mid_point - (a*(abs(n) + 1))/2), 2))
                    c1 = np.append(c1, np.around((mid_point + (a*(abs(n) + 1)/2)), 2))
            c1_p1 = np.sqrt(b_dist*b_dist + n*n*a*a)
            r1 = np.append(r1, c1_p1)
            c1_p2 = np.sqrt(b_dist*b_dist + (a*a*(n + 1)*(n + 1)))
            r2 = np.append(r2, c1_p2)
            c2 = np.append(c2,-999999)
        elif array_type in (11,15):
            c1 = np.append(c1, float(vals[1]))
            c1_z = np.append(c1_z, float(vals[2]))
            c2 = np.append(c2, float(vals[3]))
            c2_z = np.append(c2_z, float(vals[4]))
            p1 = np.append(p1, float(vals[5]))
            p1_z = np.append(p1_z, float(vals[6]))
            p2 = np.append(p2, float(vals[7]))
            p2_z = np.append(p2_z, float(vals[8]))
            pa = np.append(pa, float(vals[9]))
            elecs_all_x = np.concatenate((c1,c2,p1,p2))
            elecs_all_z = np.concatenate((c1_z,c2_z,p1_z,p2_z))
            elecs_all = np.unique(np.column_stack((elecs_all_x,elecs_all_z)), axis=0)
            if np.sum(elecs_all[:,1]) != 0: topo_flag_GA = True
            if ip_flag:
                ip = np.append(ip, float(vals[10]))
           
    #convert apparent resistivity back in to transfer resistance and vice versa
    if array_type == 2 and factor_used == False:
        K = 2*np.pi*r1
    else:
        K = geom_fac(c1, c2, p1, p2)
        
    if meas_type_flag == 'appRes': # TODO: let Survey() take care of K
        data_dict['resist'] = pa/K            
        data_dict['Rho'] = pa
    else:
        data_dict['resist'] = pa            
        data_dict['Rho'] = pa*K
        
    data_dict['dev'] = [0]*num_meas
#        data_dict['resist'].append(R)
    if ip_flag == True:
        data_dict['ip'] = ip
    else:
        data_dict['ip'] = [0]*num_meas

    if array_type == 2 and factor_used == True:
        for k in range(num_meas):
            if c2[k] > 2.5 * r1[k]:
                c2[k] = -999999
            if p2[k] > 2.5 * r1[k]:
                p2[k]= 999999

    if array_type == 8:
        total_x = np.append(total_x, c1)
        total_x = np.append(total_x, p1)
        ex_pos = np.unique(total_x)
        ey_pos = []
        for i in range(len(ex_pos)):
            ey_pos = np.append(ey_pos, b_dist)
        elec = np.column_stack((ex_pos, ey_pos))
        ex_pos = np.copy(elec)
    else:
        total_x = np.append(total_x, c1)
        total_x = np.append(total_x, c2)
        total_x = np.append(total_x, p1)
        total_x = np.append(total_x, p2)
        #convert the x electrode coordinates into indexes?
        ex_pos = np.unique(total_x)#;print(ex_pos)
    
    largo = len(c1)
    e_idx_c1 = []
    e_idx_c2 = []
    e_idx_p1 = []
    e_idx_p2 = []

    e_idx_c1 = [np.where(ex_pos == c1[i])[0][0] for i in range(largo)]
    e_idx_c1 = np.add(e_idx_c1, 1)
    data_dict['a'] = np.copy(e_idx_c1)
 
    e_idx_c2 = [np.where(ex_pos == c2[i])[0][0] for i in range(largo)]
    e_idx_c2 = np.add(e_idx_c2, 1)
    data_dict['b'] = np.copy(e_idx_c2)

    e_idx_p1 = [np.where(ex_pos == p1[i])[0][0] for i in range(largo)]
    e_idx_p1 = np.add(e_idx_p1, 1)
    data_dict['m'] = np.copy(e_idx_p1)
 
    e_idx_p2 = [np.where(ex_pos == p2[i])[0][0] for i in range(largo)]
    e_idx_p2 = np.add(e_idx_p2, 1)
    data_dict['n'] = np.copy(e_idx_p2)
    
    if array_type == 8:
        e_idx_c2_1 = [ex_pos[:,0] == c2[x] for x in range(largo)]            
        e_idx_c2_2 = [ex_pos[:, 1] == c2_y[y] for y in range(largo)]
        e_idx_c2_3 = [e_idx_c2_1[i] & e_idx_c2_2[i] for i in range(largo)]
        e_idx_c2 = [np.where(e_idx_c2_3[i])[0][0] for i in range(largo)]
        e_idx_c2 = np.add(e_idx_c2, 1)
        data_dict['b'] = np.copy(e_idx_c2)
    
        e_idx_p2_1 = [ex_pos[:,0] == p2[x] for x in range(largo)]            
        e_idx_p2_2 = [ex_pos[:, 1] == p2_y[y] for y in range(largo)]
        e_idx_p2_3 = [e_idx_p2_1[i] & e_idx_p2_2[i] for i in range(largo)]
        e_idx_p2 = [np.where(e_idx_p2_3[i])[0][0] for i in range(largo)]
        e_idx_p2 = np.add(e_idx_p2, 1)
        data_dict['n'] = np.copy(e_idx_p2)
    
    num_elec = len(ex_pos)
        
    fmt_flag = True
            
    topo_flag_idx = idx_oi + num_meas
    try:
        int(dump[topo_flag_idx])#hot fix
    except ValueError:
        topo_flag_idx+=1
    
    if int(dump[topo_flag_idx]) == 2 :#if we have topography then we should read it into the API
        topo_flag = True
        num_elec_topo =  int(dump[topo_flag_idx+1])
        ex_pos_topo=[0]*num_elec_topo
        ez_pos_topo=[0]*num_elec_topo 
        ey_pos=[0]*num_elec # actaully we can't have a y coordinate for 2d data so these will remain as zero
        ez_pos=[0]*num_elec 
        
        for i in range(num_elec_topo):
            e_pos_topo_str = dump[topo_flag_idx+2+i].strip()
            e_pos_topo_vals = re.findall(r'[-+]?\d*\.\d+|\d+', e_pos_topo_str)
#            e_pos_topo_vals = re.split(';|,|, | , | |    |\t', e_pos_topo_str)
            ex_pos_topo[i] = float(e_pos_topo_vals[0])
            ez_pos_topo[i] = float(e_pos_topo_vals[1])
            
        # finding common topography points
        elecdf = pd.DataFrame()
        elecdf['x'] = ex_pos
        elecdf['z_i'] = ez_pos
        
        elecdf_topo = pd.DataFrame()
        elecdf_topo['x'] = ex_pos_topo
        elecdf_topo['z_topo'] = ez_pos_topo
        
        if len(elecdf) != len(elecdf_topo):
            elecdf_merged = pd.merge(elecdf.copy(), elecdf_topo.copy(), how='left', on=['x'])
            ez_pos = np.array(elecdf_merged['z_topo'])
        else:
            ex_pos = ex_pos_topo.copy()
            ez_pos = ez_pos_topo.copy()        
       
    #add some protection against a dodgey file 
    if fmt_flag is False:
        raise ImportError("Error importing res2dinv input file:"+file_path+"\n the file is either unrecognised or unsupported")        
    
    #now we have indexed electrode coordinates in ex_pos :) 
    if topo_flag_GA: # if we have integrated topography in the general arrays
        ey_pos=[0]*num_elec # actaully we can't have a y coordinate for 2d data so these will remain as zero
        ez_pos=elecs_all[:,1]    
        
    elif not topo_flag: # then we dont have any topography and the electrode positions are simply given by thier x coordinates
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
    """Read in .stg file from sting (2D and 3D)
    """
    df_raw = pd.read_csv(fname, skipinitialspace=True, skiprows=3, header=None)
    elec_x =  np.concatenate(df_raw.iloc[:,[9,12,15,18]].values)
    elec_y =  np.concatenate(df_raw.iloc[:,[10,13,16,19]].values)
    elec_z =  np.concatenate(df_raw.iloc[:,[11,14,17,20]].values)
    elec_raw = np.unique(np.column_stack((elec_x,elec_y,elec_z)), axis=0)
#    elec_raw = elec_raw[elec_raw[:,1].argsort(kind='mergesort')] # to make an stable mesh, but its not working with the 3D example
    
    #detect 2D or 3D
    survey_type = '2D' if len(np.unique(elec_raw[:,1])) == 1 else '3D'
    
    if survey_type == '2D':
        elec = elec_raw[elec_raw[:,0].argsort(kind='mergesort')]# final electrode array
        a_f = [0]*len(df_raw)
        b_f = [0]*len(df_raw)
        n_f = [0]*len(df_raw)
        m_f = [0]*len(df_raw)
        
        for i in range(len(df_raw)):
            a_f[i] = list(elec[:,0]).index(df_raw.iloc[i,[9]].values)+1
            b_f[i] = list(elec[:,0]).index(df_raw.iloc[i,[12]].values)+1
            m_f[i] = list(elec[:,0]).index(df_raw.iloc[i,[15]].values)+1
            n_f[i] = list(elec[:,0]).index(df_raw.iloc[i,[18]].values)+1
    
    else: # below assumes the array is organized in a grid and not rangom XYZ values
        elecdf = pd.DataFrame(elec_raw[elec_raw[:,1].argsort(kind='mergesort')]).rename(columns={0:'x',1:'y',2:'z'})
        # organize 3D electrodes
        elecdf_groups = elecdf.groupby('y', sort=False, as_index=False)
        elecdf_lines = [elecdf_groups.get_group(x) for x in elecdf_groups.groups]
        i = 1 # index of odd lines
        while i <= len(elecdf_lines): # electrodes are laied out like a snake - although not sure if this is correct
            elecdf_lines[i]['x'] = elecdf_lines[i]['x'].values[::-1]
            i = 2*i + 1
        
        elec = np.concatenate(elecdf_lines) # final electrode array
        
        lines = np.unique(elecdf.y) # basically saying what is the y val of each line    
        
        # for final array
        a_f = []
        b_f = []
        m_f = []
        n_f = []
        
        # positions of ABMN
        array_A = df_raw.iloc[:,9:12].rename(columns={9:'x',10:'y',11:'z'})
        array_B = df_raw.iloc[:,12:15].rename(columns={12:'x',13:'y',14:'z'})
        array_M = df_raw.iloc[:,15:18].rename(columns={15:'x',16:'y',17:'z'})
        array_N = df_raw.iloc[:,18:21].rename(columns={18:'x',19:'y',20:'z'})
        
        # building A locs/labels
        array_A_groups = array_A.groupby('y', sort=False, as_index=False)
        array_A_lines = [array_A_groups.get_group(x) for x in array_A_groups.groups]
        # which lines
        for line in array_A_lines:
            line_num = np.where(lines == line['y'].iloc[0])[0][0]
            a = [0]*len(line)
            for i in range(len(line)):
                a[i] = elecdf_lines[line_num]['x'][elecdf_lines[line_num]['x'] == line['x'].iloc[i]].index[0] + 1
            a_f.extend(a)
        
        # building B locs/labels
        array_B_groups = array_B.groupby('y', sort=False, as_index=False)
        array_B_lines = [array_B_groups.get_group(x) for x in array_B_groups.groups]
        # which lines
        for line in array_B_lines:
            line_num = np.where(lines == line['y'].iloc[0])[0][0]
            b = [0]*len(line)
            for i in range(len(line)):
                b[i] = elecdf_lines[line_num]['x'][elecdf_lines[line_num]['x'] == line['x'].iloc[i]].index[0] + 1
            b_f.extend(b)
        
        # building M locs/labels
        array_M_groups = array_M.groupby('y', sort=False, as_index=False)
        array_M_lines = [array_M_groups.get_group(x) for x in array_M_groups.groups]
        # which lines
        for line in array_M_lines:
            line_num = np.where(lines == line['y'].iloc[0])[0][0]
            m = [0]*len(line)
            for i in range(len(line)):
                m[i] = elecdf_lines[line_num]['x'][elecdf_lines[line_num]['x'] == line['x'].iloc[i]].index[0] + 1
            m_f.extend(m)
        
        # building N locs/labels
        array_N_groups = array_N.groupby('y', sort=False, as_index=False)
        array_N_lines = [array_N_groups.get_group(x) for x in array_N_groups.groups]
        # which lines
        for line in array_N_lines:
            line_num = np.where(lines == line['y'].iloc[0])[0][0]
            n = [0]*len(line)
            for i in range(len(line)):
                n[i] = elecdf_lines[line_num]['x'][elecdf_lines[line_num]['x'] == line['x'].iloc[i]].index[0] + 1
            n_f.extend(n)
    
    #build df
    df = pd.DataFrame()
    df['a'] = np.array(a_f)
    df['b'] = np.array(b_f)
    df['n'] = np.array(n_f)
    df['m'] = np.array(m_f)
    df['resist']=df_raw.iloc[:,4]

    #detecting IP (col 21 would be floats of IP vals)
    try:
        float(df_raw.iloc[0,21])
        df['ip']=df_raw.iloc[0,21]
    except:
        df['ip']=[0]*len(df_raw)
    
    #for pole-pole and pole-dipole arrays
    elec[elec > 9999] = 999999
    elec[elec < -9999] = -999999
    df = df.query('a!=b & b!=m & m!=n & a!=m & a!=n & b!=n') # removing data where ABMN overlap
    
    return elec,df


#%% ABEM Lund parsers (2D and 3D)
def ericParser(file_path):
    """
    Reads *.ohm ASCII-files with information related to the profile, comment,
    first station coordinate, date and time, version of data collection program
    used, electrode take-out spacing and the protocol files used.
    3D *.ohm files should use array type = 12 (tomography). type "b dist="num
    in the first row of the *.ohm file. If no b dist exist then it is assumed
    that b dist= distance among electrodes. There num = distance among parallel
    lines of electrodes. 
    """
    #ericParser Rev 2020-05-18
    try: # see if it's 3D
        fh = open(file_path,'r')#open file handle for reading
        dump = fh.readlines()#cache file contents into a list
        fh.close()#close file handle, free up resources
        
        #declaration of variables    
        proto_file = []   
        array_type = []
        proto_org = []
        num_meas = []
        mid_st_coord = []
        idx_meas = []
        c1 = np.array(())
        c2 = np.array(())
        p1 = np.array(())
        p2 = np.array(())
        pt = np.array(())
        pa = np.array(())
        var_coef = np.array(())
        n_cycles = np.array(())
        n_tot_cycles = np.array(())
        total_x = np.array(())
        h_dist = np.array(())
        a_dist = np.array(())
        b_dist = np.array(())
        e_x = np.array(())
        e_y = np.array(())
        e_z = np.array(())    
    #    data_dict = {'a':[],'b':[],'m':[],'n':[],'Rho':[],'ip':[],'resist':[],'dev':[]}
        df = pd.DataFrame()
        
        tot_num_meas = 0
        #first find the general information
        idx_oi = 0
        line = dump[idx_oi]
        sur_name = line.strip() #name of survey
        if "b dist=" in sur_name:
            vals = line.strip().split()
            b_dist = float(vals[2])
        idx_oi += 1
        line = dump[idx_oi]
        vals = line.strip().split()
        x_location = float(vals[0])#First midstation coordinate
        idx_oi += 1
        line = dump[idx_oi]
        vals = line.strip().split()
        date_time_sur = str(vals[0]) + str('  ') + str(vals[1]) 
        eric_version = str(vals[2]) + str(': ') + str(vals[3])
        idx_oi += 1
        line = dump[idx_oi]
        vals = line.strip().split()
        a_spac = float(vals[0]) #electrode spacing
        idx_oi += 1
        no_protocols = 0
        idx_proto_file = 1
        idx_oi += 1
        line = dump[idx_oi]
        vals = line.strip().split()
        elec_cable = int(vals[0][2:4])    
    #    no_cables = int(vals[0][-2:-1])    
        no_cables = int([s for s in vals[0] if s.isdigit()].pop())
        proto_file.append(str(vals[0]))
        array_type.append(int(vals[1]))
        pdip_pp_flag = False
        polpol_flag = False
            
        #First find how many protocol, measurements and mid first location are 
        #included in the *.OHM file
        
        for i, line in enumerate(dump):
            pro_type = line.strip().split('.')
            if 'ORG' in pro_type or 'UP' in pro_type or 'DWN' in pro_type:
                proto_org.append(str(line)) 
                no_protocols = no_protocols + 1
                linea = dump[i+1]
                vals = linea.strip().split()
                num_meas.append(int(vals[0]))
                mid_st_coord.append(float(dump[i+2]))
                idx_oi = i + 4
                idx_meas.append(idx_oi)
                
        for i in range(len(num_meas)):
            for k in range(num_meas[i]):
                line = dump[idx_meas[i] + k]
                vals = line.strip().split()
                c1 = np.append(c1, float(vals[0]))
                c2 = np.append(c2, float(vals[1]))
                p1 = np.append(p1, float(vals[2]))
                p2 = np.append(p2, float(vals[3]))
                pt = np.append(pt, float(vals[4]))
                var_coef = np.append(var_coef, float(vals[5]))
                n_cycles = np.append(n_cycles, int(vals[6]))
                n_tot_cycles = np.append(n_tot_cycles, int(vals[7]))
    
        if array_type == [12]:
            if sur_name == "":
                b_dist = a_spac
            e_x = np.linspace(0, (no_cables - 1)*b_dist, num= no_cables)
            e_y = np.linspace(0, (elec_cable - 1)*a_spac, num= elec_cable)
            e_z = [0]
            elec = np.vstack((ndmesh(e_x,e_y,e_z))).reshape(3,-1).T
            num_elec = len(elec)
            df['a'] = np.copy(c1)
            df['b'] = np.copy(c2)
            df['m'] = np.copy(p1)
            df['n'] = np.copy(p2)
            k = geom_factor_3D(df, elec, array_type)
        else:
            min_dist_c1 = min(c1)
            min_dist_p1 = min(p1)
            if min_dist_c1 <= min_dist_p1:
                min_dist = min_dist_c1
            else:
                min_dist = min_dist_p1
                
            max_dist_c2 = max(c2)
            max_dist_p1 = max(p1)
               
            if max_dist_c2 >= max_dist_p1:
                max_dist = max_dist_c2
            else:
                max_dist = max_dist_p1
            
            if min_dist <= 0.0:
                half_dist = abs(min_dist)
            else:
                half_dist = 0.0
                
            max_dist_p2 = max(p2) 
            largo = len(c1)    
            for k in range(largo):
                h_dist = np.append(h_dist, half_dist)
            
            c1 = np.add(c1, h_dist)
            if max_dist_c2 == 1e+38:
                for k in range(largo):
                    c2[k] = -999999
            else:
                c2 = np.add(c2, h_dist)
            p1 = np.add(p1, h_dist)
            if max_dist_p2 == 1e+38:
                for k in range(largo):
                    p2[k] = 999999
            else:
                p2 = np.add(p2, h_dist)
            
            total_x = np.append(total_x, c1)
            total_x = np.append(total_x, c2)
            total_x = np.append(total_x, p1)
            total_x = np.append(total_x, p2)
            ex_pos = np.unique(total_x)
            
            num_elec = len(ex_pos)
            e_idx_c1 = []
            e_idx_c2 = []
            e_idx_p1 = []
            e_idx_p2 = []
        
            e_idx_c1 = [np.where(ex_pos == c1[i])[0][0] for i in range(largo)]
            e_idx_c1 = np.add(e_idx_c1, 1)
            df['a'] = np.copy(e_idx_c1)
        
            e_idx_c2 = [np.where(ex_pos == c2[i])[0][0] for i in range(largo)]
            e_idx_c2 = np.add(e_idx_c2, 1)
            df['b'] = np.copy(e_idx_c2)
            e_idx_p1 = [np.where(ex_pos == p1[i])[0][0] for i in range(largo)]
            e_idx_p1 = np.add(e_idx_p1, 1)
            df['m'] = np.copy(e_idx_p1)
            e_idx_p2 = [np.where(ex_pos == p2[i])[0][0] for i in range(largo)]
            e_idx_p2 = np.add(e_idx_p2, 1)   
            df['n'] = np.copy(e_idx_p2)
            df['resist'] = np.copy(pt)
                    
            k = geom_fac(c1, c2, p1, p2)
            ey_pos=[0]*num_elec
            ez_pos=[0]*num_elec  
            elec = np.column_stack((ex_pos,ey_pos,ez_pos))
            
            #for pole-pole and pole-dipole arrays
            elec[elec > 9999] = 999999
            elec[elec < -9999] = -999999
    
       
        df['resist'] = np.copy(pt)
        df['Rho'] = pt*k 
        df['dev'] = (var_coef * n_tot_cycles * pt)/100
        df['ip'] = [0]*len(c1)
        #we dont have any topography at x coordinates
        array = df[['a','b','m','n']].values
        arrayMin = np.min(np.unique(np.sort(array.flatten())))
        if arrayMin != 0: # all surveys must start from x = 0
            array -= arrayMin
        df[['a','b','m','n']] = (array+1).astype(int)
           
        df = df[['a','b','m','n','Rho','dev','ip','resist']] # reorder columns to be consistent with the syscal parser
    except:
        elec, df = ericParser2D(file_path) # well it seems to be 2D
        
    return elec,df

def ericParser2D(file_path):
    """
    Reads *.ohm ASCII-files with information related to the profile, comment,
    first station coordinate, date and time, version of data collection program
    used, electrode take-out spaing and the protocol files used.
    """
    
    fh = open(file_path,'r')#open file handle for reading
    dump = fh.readlines()#cache file contents into a list
    fh.close()#close file handle, free up resources
    
    #declaration of variables    
    proto_file = []   
    array_type = []
    proto_org = []
    num_meas = []
    mid_st_coord = []
    idx_meas = []
    c1 = np.array(())
    c2 = np.array(())
    p1 = np.array(())
    p2 = np.array(())
    pt = np.array(())
    pa = np.array(())
    var_coef = np.array(())
    n_cycles = np.array(())
    n_tot_cycles = np.array(())
    total_x = np.array(())
    h_dist = np.array(())
    data_dict = {'a':[],'b':[],'m':[],'n':[],'Rho':[],'ip':[],'resist':[],'dev':[]}
    tot_num_meas = 0
    #first find the general information
    idx_oi = 0
    line = dump[idx_oi]
    sur_name = line.strip() #name of survey
    idx_oi += 1
    line = dump[idx_oi]
    vals = line.strip().split()
    x_location = float(vals[0])#First midstation coordinate
    idx_oi += 1
    line = dump[idx_oi]
    vals = line.strip().split()
    date_time_sur = str(vals[0]) + str('  ') + str(vals[1]) 
    eric_version = str(vals[2]) + str(': ') + str(vals[3])
    idx_oi += 1
    line = dump[idx_oi]
    vals = line.strip().split()
    a_spac = float(vals[0]) #electrode spacing
    idx_oi += 1
    #line = dump[idx_oi]
    #vals = line.strip().split()
    #no_protocols = int(vals[0])#no. of protocol used
    no_protocols = 0
    idx_proto_file = 1
    idx_oi += 1
    line = dump[idx_oi]
    vals = line.strip().split()
    proto_file.append(str(vals[0]))
    array_type.append(int(vals[1]))
        
    #First find how many protocol, measurements and mid first location are 
    #included in the *.OHM file
    
    for i, line in enumerate(dump):
        pro_type = line.strip().split('.')
        if 'ORG' in pro_type or 'UP' in pro_type or 'DWN' in pro_type:
            proto_org.append(str(line)) 
            no_protocols = no_protocols + 1
            linea = dump[i+1]
            vals = linea.strip().split()
            num_meas.append(int(vals[0]))
            mid_st_coord.append(float(dump[i+2]))
            idx_oi = i + 4
            idx_meas.append(idx_oi)
            
    for i in range(len(num_meas)):
        for k in range(num_meas[i]):
            line = dump[idx_meas[i] + k]
            vals = line.strip().split()
            c1 = np.append(c1, float(vals[0]))
            c2 = np.append(c2, float(vals[1]))
            p1 = np.append(p1, float(vals[2]))
            p2 = np.append(p2, float(vals[3]))
            pt = np.append(pt, float(vals[4]))
            var_coef = np.append(var_coef, float(vals[5]))
            n_cycles = np.append(n_cycles, int(vals[6]))
            n_tot_cycles = np.append(n_tot_cycles, int(vals[7]))
            data_dict['ip'].append(0)
            
    min_dist_c1 = min(c1)
    min_dist_p1 = min(p1)
    if min_dist_c1 <= min_dist_p1:
        min_dist = min_dist_c1
    else:
        min_dist = min_dist_p1
        
    max_dist_c2 = max(c2)
    max_dist_p1 = max(p1)
       
    if max_dist_c2 >= max_dist_p1:
        max_dist = max_dist_c2
    else:
        max_dist = max_dist_p1
    
    if min_dist <= 0.0:
        half_dist = abs(min_dist)
    else:
        half_dist = 0.0
        
    max_dist_p2 = max(p2) 
    largo = len(c1)    
    for k in range(largo):
        h_dist = np.append(h_dist, half_dist)
    
    c1 = np.add(c1, h_dist)
    if max_dist_c2 == 1e+38:
        for k in range(largo):
            c2[k] = -999999
    else:
        c2 = np.add(c2, h_dist)
    p1 = np.add(p1, h_dist)
    if max_dist_p2 == 1e+38:
        for k in range(largo):
            p2[k] = 999999
    else:
        p2 = np.add(p2, h_dist)
    
    total_x = np.append(total_x, c1)
    total_x = np.append(total_x, c2)
    total_x = np.append(total_x, p1)
    total_x = np.append(total_x, p2)
    ex_pos = np.unique(total_x)
    
    num_elec = len(ex_pos)
    e_idx_c1 = []
    e_idx_c2 = []
    e_idx_p1 = []
    e_idx_p2 = []
    e_idx_c1 = [np.where(ex_pos == c1[i])[0][0] for i in range(largo)]
    e_idx_c2 = [np.where(ex_pos == c2[i])[0][0] for i in range(largo)]
    e_idx_p1 = [np.where(ex_pos == p1[i])[0][0] for i in range(largo)]
    e_idx_p2 = [np.where(ex_pos == p2[i])[0][0] for i in range(largo)]
    
    e_idx_c1 = np.add(e_idx_c1, 1)
    e_idx_c2 = np.add(e_idx_c2, 1)
    e_idx_p1 = np.add(e_idx_p1, 1)
    e_idx_p2 = np.add(e_idx_p2, 1)
    
    data_dict['a'] = np.copy(e_idx_c1)
    data_dict['b'] = np.copy(e_idx_c2)
    data_dict['m'] = np.copy(e_idx_p1)
    data_dict['n'] = np.copy(e_idx_p2)
    data_dict['resist'] = np.copy(pt)
            
    AM = np.abs(c1-p1)
    BM = np.abs(c2-p1)
    AN = np.abs(c1-p2)
    BN = np.abs(c2-p2)
    K = 2*np.pi/((1/AM)-(1/BM)-(1/AN)+(1/BN))#geometric factor
    data_dict['Rho'] = pt*K
    data_dict['dev'] = (var_coef * n_tot_cycles * pt)/100
    #we dont have any topography at x coordinates
    ey_pos=[0]*num_elec
    ez_pos=[0]*num_elec  
    elec = np.column_stack((ex_pos,ey_pos,ez_pos))
    
    #for pole-pole and pole-dipole arrays
    elec[elec > 9999] = 999999
    elec[elec < -9999] = -999999
       
    df = pd.DataFrame(data=data_dict) # make a data frame from dictionary
    df = df[['a','b','m','n','Rho','dev','ip','resist']] # reorder columns to be consistent with the syscal parser
    
    return elec, df


#%% 
def lippmannParser(fname):
    """Read in *.tx0 file from Lippmann instruments
    """
    with open(fname, 'r') as fh:
        dump = fh.readlines()

    #getting electrode locations
    elec_lineNum_s = [i for i in range(len(dump)) if '* Electrode positions *' in dump[i]]
    elec_lineNum_e = [i-1 for i in range(len(dump)) if '* Remote electrode positions *' in dump[i]]
    elec_nrows = elec_lineNum_e[0] - elec_lineNum_s[0]
    elec_raw = pd.read_csv(fname, delim_whitespace=True, skiprows=elec_lineNum_s[0]+1, nrows=elec_nrows, header=None)
    elec = np.array(elec_raw.iloc[:,-3:])
    
    #for pole-pole and pole-dipole arrays
    elec[elec > 9999] = 999999
    elec[elec < -9999] = -999999

    #getting data
    data_linNum_s = [i for i in range(len(dump)) if '* Data *********' in dump[i]]
    data_headers = dump[data_linNum_s[0]+1].split()[1:]
    df = pd.read_csv(fname, delim_whitespace=True, skiprows=data_linNum_s[0]+3, names=data_headers).drop('n', axis=1) # don't know what this "n" is!!
    df = df.rename(columns={'A':'a',
                            'B':'b',
                            'M':'m',
                            'N':'n',
                            'I':'i',
                            'U':'vp'})
    if 'phi' in df.columns:
        df = df.rename(columns={'phi':'ip'})
        df = df[['a','b','m','n','i','vp','ip']]
    else:
        df = df[['a','b','m','n','i','vp']]
        df['ip'] = 0
        
    df = df.query("i != '-' & vp != '-' & ip != '-'").astype(float)    
    
    #calculations
    df['resist'] = df['vp']/df['i']
    df = df.astype({'a':int, 'b':int, 'm':int, 'n':int})
    
    return elec, df

# elec, df = lippmannParser(testdir + 'parser/Lippmann_1.tx0')

#%%
def aresParser(fname, spacing=None):
    """Read in *.2dm file from ARES II
    """
    with open(fname, 'r') as fh:
        dump = fh.readlines()
    
    #getting spacing
    spacing_lineNum = [i for i in range(len(dump)) if 'Electrode distance' in dump[i]]
    if spacing_lineNum != []:
        spacing = dump[spacing_lineNum[0]].split()[2]

    #getting data
    data_linNum_s = [i for i in range(len(dump)) if 'Measured data' in dump[i]]
    df = pd.read_csv(fname, delim_whitespace=True, skiprows=data_linNum_s[0]+1, index_col=False)
    df = df.rename(columns={'C1[el]':'a',
                            'C2[el]':'b',
                            'P1[el]':'m',
                            'P2[el]':'n',
                            'I[mA]':'i',
                            'U[mV]':'vp'})
    
    #Building electrode locations
    df[['a','b','m','n']] = df[['a','b','m','n']].apply(pd.to_numeric, errors='coerce') # there are strange eletrode conbinations sometimes (e.g., 14*1)
    df = df.dropna()
    array = df[['a','b','m','n']].values
    
    # arrayMin = np.min(np.unique(np.sort(array.flatten())))
    # if arrayMin != 0: # all surveys must start from x = 0
    #     array -= arrayMin
    val = np.sort(np.unique(array.flatten())) # unique electrodes positions/labels required are 
    # elecLabel = 1 + np.arange(len(val))
    # newval = elecLabel[np.searchsorted(val, array)] # magic ! https://stackoverflow.com/questions/47171356/replace-values-in-numpy-array-based-on-dictionary-and-avoid-overlap-between-new
    # df.loc[:,['a','b','m','n']] = newval
    
    #calculations
#    if 'EP[mV]' in df.columns: # TODO: correct this when you figure out how the IP values are calculated EP may acctualy be SP!!
#        df['ip'] = df['EP[mV]']/df['vp'] # not sure this is correct way of calculating IP
#        df = df[['a','b','m','n','i','vp','ip']]
#    else:
    
    
    # finding IP columns
    
    ip_cols = [col for col in df.columns if all(i in col for i in ['IP', '[%]'])] # assuming IP columns are those with IP and [%]
    if ip_cols!= []:
        # df[ip_cols] = df[ip_cols]/100 # not sure about this... who reports IP in % anyway?!
        df['ip'] = df[ip_cols].mean(axis=1) # average of all IP windows
        df = df[['a','b','m','n','i','vp','ip']]
        
    else:
        df = df[['a','b','m','n','i','vp']] # should be under "else:"
        df['ip'] = 0 # should be under "else:"
    
    df = df.query("i != '-' & vp != '-' & ip != '-'").astype(float)    
    
    df['resist'] = df['vp']/df['i']
    
    #builting electrode table
    if spacing is not None:
        elec = np.c_[val*float(spacing), np.zeros((len(val),2))]
    else:
        elec = np.c_[val, np.zeros((len(val),2))]
    
    #for pole-pole and pole-dipole arrays
    elec[elec > 9999] = 999999
    elec[elec < -9999] = -999999
    
    return elec, df


#%% BERT format parser  
def bertParser(fname):
    f = open(fname, "r")
    
    dump = f.readlines()
    line = 0
    
    # skip comment lines
    while dump[line][0] == '#':
        line += 1

    numElec = re.findall(r'[-+]?\d*\.\d+|\d+', dump[line])
    if len(numElec) == 1: # we have number of elecs
        line += 1
    
    elecHeaders = re.findall(r'#|x|y|z', dump[line])
    if len(elecHeaders) != 0: # we have elec location headers (x, y, z)
        line += 1
    
    elec_list = []
    elecLocs0 = re.findall(r'[-+]?\d*\.\d+|\d+', dump[line])
    elecLocs_line = elecLocs0.copy()
    while len(elecLocs_line) == len(elecLocs0):
        elecLocs_line = re.findall(r'[-+]?\d*\.\d+|\d+', dump[line])
        elec_list.append(elecLocs_line)
        line += 1
    
    elec = np.array(elec_list[:-1]).astype(float)
    
    if elec.shape[1] != 3: # we have xz format so conver into xyz
        elec = np.c_[elec[:,0], np.zeros(len(elec)), elec[:,1]]
    
    vals = re.findall(r'[-+]?\d*\.\d+|\d+', dump[line])
    while len(vals) < 4: # finding the data line
        line += 1
        vals = re.findall(r'[-+]?\d*\.\d+|\d+', dump[line])
    
    headers = re.findall(r'[A-Za-z]+', dump[line-1]) # for finding data types

    topo_check_vals = len(re.findall(r'[-+]?\d*\.\d+|\d+', dump[line])) # TODO: is topography included without any flags?
    df_list = []
    for val in dump[line:]: # reding data
        df_list.append(re.findall(r'[-+]?\d*\.\d+|\d+', val))
    
    df = pd.DataFrame(np.array(df_list).astype(float)) # getting the electrode array
    df.columns = headers
    
    resristance_list = ['r', 'R', 'rho', 'Rho', 'RHO']    
    resistivity_list = ['rhoa', 'Rhoa', 'Ra', 'ra', 'RA', 'RHOA']    
    ip_list = ['ip', 'IP']    
    iv_list = ['i', 'I', 'u', 'U', 'i/mA', 'u/mV']    
    resErr_list = ['err', 'ERR', 'Err', 'err/%']    
    ipErr_list = ['ipErr', 'IPERR', 'IPerr', 'iperr']
    
    #check whether resistance is given or app. res or (I and V)
    if any(header in headers for header in resristance_list):
        header = [col for col in resristance_list if col in headers][0]
        df = df.rename(columns = {header:'resist'})

    elif any(header in headers for header in resistivity_list):
        header = [col for col in resistivity_list if col in headers][0]
        df = df.rename(columns = {header:'app'}) # no resistance calc required here, Survey() takes care of it

    elif any(header in headers for header in iv_list):
        i = [col for col in ['i', 'I', 'i/mA'] if col in headers][0]
        v = [col for col in ['u', 'U', 'u/mV'] if col in headers][0]
        df = df.rename(columns = {i:'i', v:'vp'})
        df['resist'] = df['vp']/df['i']

    else:
        raise ValueError('Data does not contain enough columns for obtaining resistance/resistivity')
        
    if any(header in headers for header in ip_list): # IP check
        header = [col for col in ip_list if col in headers][0]
        df = df.rename(columns = {header:'ip'})
        
    if any(header in headers for header in resErr_list): # R_err check (is it stacking error or magnitude error?)
        header = [col for col in resErr_list if col in headers][0]
        if '%' in header: # error should be in fractions
            df[header].values /= 100
        df = df.rename(columns = {header:'magErr'})
        
    if any(header in headers for header in ipErr_list): # IP_err check
        header = [col for col in ipErr_list if col in headers][0]
        df = df.rename(columns = {header:'phiErr'})
    
    if 'ip' not in df.columns:
        df['ip'] = np.nan
    df = df.rename(columns = dict(zip(df.columns[:4], ['a', 'b', 'm', 'n']))) # make sure headers are a b m n
    df = df.astype({'a':int, 'b':int, 'm':int, 'n':int})
    
    return elec,df

#%% E4D .srv format 
def srvParser(fname):
    fh = open(fname,'r')
    # electrodes come first in the file 
    nelec = int(fh.readline().split()[0]) # number of electrode coordinates 
    # now read in electrode coords
    elec = np.zeros((nelec,3))
    for i in range(nelec):
        line = fh.readline().split()
        elec[i,0] = float(line[1])
        elec[i,1] = float(line[2])
        elec[i,2] = float(line[3])
    #read in measurements 
    line = fh.readline()
    while len(line.split())==0:
        line = fh.readline()
        
    nmeas = int(line.split()[0])
    a = [0]*nmeas
    b = [0]*nmeas
    n = [0]*nmeas
    m = [0]*nmeas
    Tr = [0]*nmeas
    err = [0]*nmeas
    for i in range(nmeas):
        line = fh.readline().split()
        a[i] = int(line[1])
        b[i] = int(line[2])
        m[i] = int(line[3])
        n[i] = int(line[4])
        Tr[i] = float(line[5])
        if len(line) == 7:
            err[i] = float(line[6]) # grab the error model from the file? 
            
    if any(err) != 0:
        print('E4D error column found')
    else:
        err = [np.nan]*nmeas
        
    #put data into correct format
    data_dict = {}
    data_dict['a']=a
    data_dict['b']=b
    data_dict['n']=n
    data_dict['m']=m
    data_dict['resist']=Tr
    data_dict['magErr']=err
    data_dict['Rho']=[0]*nmeas
    data_dict['dev']=[0]*nmeas
    data_dict['ip']=[0]*nmeas
    df = pd.DataFrame(data=data_dict) # make a data frame from dictionary
    df = df[['a','b','m','n','Rho','dev','ip','resist','magErr']] # reorder columns to be consistent with the syscal parser
    
    return elec, df 
