# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 11:12:08 2019

@author: Sina
"""
import os, io 
import numpy as np
import pandas as pd
pd_vns = pd.__version__ # string version 
pd_vn = float(pd_vns.split('.')[0]) + float('0.'+pd_vns.split('.')[1]) # version as float 

def write2Res2DInv(param, fname, df, elec, typ='R2'):
    """Writes a Res2DInv format file.
    
    Parameters
    ----------
    param : dict
        Dictionnary of parameters to be used.
    fname : str
        Path of the file to be saved.
    df : DataFrame
        DataFrame containing measurements
    elec : Array
        Array containing topography information
    typ : str
        Type of file either `R2`, `cR2`, `R3t`, `cR3t`.
    
    Returns
    -------
    String to be writting in the ".dat".
    """
    dparam = {
            'lineTitle':'My beautiful survey',
            'spacing':1.0,
            'array_type':11,
            'array_spec':0,
            'header':'Type of measurement (0=app. resistivity,1=resistance)',
            'res_type':1, # default resistance
            'num_meas':0, #number of measurements - to be determined
            'x_loc_type':2, # != 0 default, 2 for location when there is topography
            'ip_flag':0, # 0 no IP, 1 with IP.
            'ip_type':'Chargeability', # IP type - default is Chargeability
            'ip_unit':'mV/V', # IP unit - default is mV/V
            'ip_spec':'0.12,1.0', # IP delay and integration time - Throw some default values!!
            'topo_header':'Topography in separate list',
            'num_topo':0, # number of topo data - to be determined
            'topo_elec_num':1, # number with first electrode
            'end_zeros':'0\n0\n0\n0' # zeros at end of file
            }
    
    # check if values are missing
    for a in dparam:
        if a not in param: # parameter missing
            param[a] = dparam[a]
            
    # create header text for .dat file
    content = ''
    content = content + '{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n'.format(
                        param['lineTitle'],
                        param['spacing'],
                        int(param['array_type']),
                        int(param['array_spec']),
                        param['header'],
                        int(param['res_type']),
                        int(param['num_meas']),
                        int(param['x_loc_type']))
    if typ == 'R2':
        param['ip_flag'] = 0
        content = content + '{}\n'.format(int(param['ip_flag']))
        
    elif typ == 'cR2':
        param['ip_flag'] = 1
        content = content + '{}\n{}\n{}\n{}\n'.format(
                            int(param['ip_flag']),
                            param['ip_type'],
                            param['ip_unit'],
                            param['ip_spec'])
        
    # formattin measurement points and adding topography info
    a = df.a#-df.a[0] # -1 for positional issues! first position should be at "0" Updat: NOT NEEDED anymore
    az = 0*a.rename('az')
    b = df.b#-df.a[0]
    bz = 0*b.rename('bz')
    m = df.m#-df.a[0]
    mz = 0*m.rename('mz')
    n = df.n#-df.a[0]
    nz = 0*n.rename('nz')
    if 'resist' in df.columns:
        resist = df.resist
    else:
        resist = df.res
    if 'ip' in df.columns:
        ip = df.ip
    elif'phase' in df.columns:
        ip = df.phase
    else:
        ip = 0*resist.rename('ip')
    dfr2d = pd.concat((a,az,b,bz,m,mz,n,nz,resist,ip), axis =1)
    param['num_meas'] = len(dfr2d)
    if typ == 'R2':
        dfr2d = dfr2d.drop(['ip'], axis=1)
        content = content + ''.join(['4 {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {}\n']*len(dfr2d)).format(
            *np.array(dfr2d).flatten())
        
    elif typ == 'cR2':
        content = content + ''.join(['4 {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {} {}\n']*len(dfr2d)).format(
            *np.array(dfr2d).flatten())
    
    if np.sum(elec[:,2]) != 0: # if we have topography
        topodf = pd.DataFrame(elec[:,[0,2]])
        param['num_topo'] = len(topodf)
        content = content + '{}\n{}\n{}\n'.format(
                            param['topo_header'],
                            int(param['x_loc_type']),
                            int(param['num_topo']))
        content = content + ''.join(['{} {}\n']*len(topodf)).format(
            *np.array(topodf).flatten())
        content = content + '{}\n'.format(int(param['topo_elec_num']))
        
    content = content + '{}'.format(param['end_zeros'])
    
#    fname = param['lineTitle']+'.dat'
    with open(fname,'w') as f:
        f.write(content)
    
    return content # if needed!


def write2csv(fname, dfi, elec, typ='R2'):
    """Writes a clean csv format file.
    
    Parameters
    ----------
    fname : str
        Path of the file to be saved.
    dfi : DataFrame
        DataFrame containing measurements
    elec : Array
        Array containing topography information
    typ : str
        Type of file either `R2`, `cR2`, `R3t`, `cR3t`.
    """
    df = dfi[['a','b','m','n','resist','ip']]
    df = df.rename(columns = {'resist':'Resistance',
                              'ip':'Chargeability'})
    if 'recipMean' in dfi.columns:
        df['Mean_R_Error'] = dfi['recipMean']
        df['Relative_R_Error'] = dfi['reciprocalErrRel']
    if 'reci_IP_err' in dfi.columns:
        df['Recipraocal_IP_Error'] = dfi['reci_IP_err']
    
    if typ  == 'R2':
        df = df.drop(['Chargeability'], axis=1)
        if 'Recipraocal_IP_Error' in df.columns:
            df = df.drop(['Recipraocal_IP_Error'], axis=1)
            
    df.to_csv(fname, index=False)
    if np.sum(elec[:,2]) != 0: # if we have topography
        topodf = pd.DataFrame(elec[:,[0,2]]).rename(columns = {0:'X [m]',1:'Z [m]'})
        topofname = fname[:-4]+'_topography'+fname[-4:]
        topodf.to_csv(topofname, index=False)
        

def writeSrv(fname, df, elec): # pragma: no cover
    """Export .srv format for which is compatible with E4D. The e4d survey
    file includes the electrode locations, in addition to the scheduling 
    matrix. 
    
    Paramters
    ------------
    fname: string, optional
        Where the output file will be written to. 
    """
    if not isinstance(fname,str):
        raise ValueError('fname must be a string')
    
    if fname is None: # rename output file name to that of the survey name
        fname = 'protocol' + '.srv'
    
    fh = open(fname,'w')
    numelec = elec.shape[0] # number of electrodes 
    fh.write('%i number of electrodes\n'%numelec)
    for i in range(numelec):
        line = '{:d} {:f} {:f} {:f} {:d}\n'.format(i+1,
                elec[i,0],#x coordinate
                elec[i,1],#y coordinate
                elec[i,2],#z coordinate
                1)#buried flag 
        fh.write(line)
    #now write the scheduling matrix to file 
    # ie = df['irecip'].values >= 0 # reciprocal + non-paired
    # df = df[ie]
    nomeas = len(df) # number of measurements 
    df = df.reset_index().copy()
    fh.write('\n%i number of measurements \n'%nomeas)
    
    if not 'resError' in df.columns or all(np.isnan(df['resError'])) is True: # the columns don't exist or are empty 
        #estimate error if not given 
        if 'resist' in df.columns:
            res = np.array(df['resist'])
        else:
            res = np.array(df['res'])
        a_wgt = 0.1
        b_wgt = 0.2
        var_res = (a_wgt*a_wgt)+(b_wgt*b_wgt) * (res*res)
        std_res = np.sqrt(var_res)
        df['resError'] = std_res
        
    
    # format >>> m_indx a b m n V/I stdev_V/I
    for i in range(nomeas): 
        line = '{:d} {:d} {:d} {:d} {:d} {:f} {:f}\n'.format(i+1,
                int(df['a'][i]),
                int(df['b'][i]),
                int(df['m'][i]),
                int(df['n'][i]),
                df['recipMean'][i],
                df['resError'][i])
        fh.write(line)
    
def to_csv(df, fh, sep=',', index=True, header=True, float_format = None, 
           lineterminator=None):
    """
    Replacement for pandas to_csv function, that doesnt throw lineterminator 
    error. 

    Parameters
    ----------
    df : Pandas.Dataframe or Dict
        DESCRIPTION.
    fh : str or TextIOWrapper
        DESCRIPTION.
    index : bool, optional
        DESCRIPTION. The default is True.
    header : bool, optional
        DESCRIPTION. The default is True. 
    lineterminator : str, optional
        DESCRIPTION. The default is os.linesep. 

    Returns
    -------
    None.

    """
    # check pandas version, if greater than 1.5.0 then use to_csv as its quicker 
    if pd_vn >= 1.5 and isinstance(df,pd.DataFrame):
        df.to_csv(fh,sep = sep, index = index, header=header, 
                  float_format = float_format, lineterminator = lineterminator)
        return 
    
    #otherwise use homebuilt solution (also works for dictionarys)
    openedfile = False 
    if isinstance(fh, str):
        fh = open(fh,'w') # then assume the string is a file path and hence
        # open a file 
        openedfile = True 
    if not isinstance(fh, io.TextIOWrapper):
        raise TypeError('Expected a TextIOWrapper object for the file handle ')
        
    if lineterminator is None:
        lineterminator = os.linesep 
        
    columns = df.keys() # works with both a dictionary and pandas data frame 
    
    ncol = len(columns) # number of columns 
    hlist = [h for h in columns]# header list 
    nline = len(df[hlist[0]])# number of lines
    
    if header:
        line = ''
        if index: 
            line += 'index' + sep 
            
        for i,h in enumerate(hlist):
            line += h
            if i == ncol-1:
                line += h + lineterminator
                break 
            line += sep 
        fh.write(line)
            
    for i in range(nline):
        line = ''
        if index:
            line += '%i'%i + sep 
        for j,h in enumerate(hlist):
            x = df[h][i]
            line += '{:}'.format(x)
            if j == ncol-1: 
                line += lineterminator 
                break
            line += sep 
        fh.write(line)
    
    if openedfile:
        fh.close() 
