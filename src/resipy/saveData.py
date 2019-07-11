# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 11:12:08 2019

@author: Sina
"""

import numpy as np
import pandas as pd
import os

def write2Res2DInv(param, dirname, df, elec, typ  = 'R2'):
    '''Writes a Res2DInv format file'''
    
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
    a = df.a-1 # -1 for positional issues! first position should be at "0"
    az = 0*a.rename('az')
    b = df.b-1
    bz = 0*b.rename('bz')
    m = df.m-1
    mz = 0*m.rename('mz')
    n = df.n-1
    nz = 0*n.rename('nz')
    resist = df.resist
    ip = df.ip
    dfr2d = pd.concat((a,az,b,bz,m,mz,n,nz,resist,ip), axis =1)
    param['num_meas'] = len(dfr2d)
    if typ == 'R2':
        dfr2d = dfr2d.drop(['ip'], axis=1)
        content = content + ''.join(['4 {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f}\n']*len(dfr2d)).format(
            *np.array(dfr2d).flatten())
        
    elif typ == 'cR2':
        content = content + ''.join(['4 {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f}\n']*len(dfr2d)).format(
            *np.array(dfr2d).flatten())
    
    if np.sum(elec[:,2]) != 0: # if we have topography
        topodf = pd.DataFrame(elec[:,[0,2]])
        param['num_topo'] = len(topodf)
        content = content + '{}\n{}\n{}\n'.format(
                            param['topo_header'],
                            int(param['x_loc_type']),
                            int(param['num_topo']))
        content = content + ''.join(['{:.2f} {:.2f}\n']*len(topodf)).format(
            *np.array(topodf).flatten())
        content = content + '{}\n'.format(int(param['topo_elec_num']))
        
    content = content + '{}'.format(param['end_zeros'])
    
    fname = param['lineTitle']+'.dat'
    with open(os.path.join(dirname,fname),'w') as f:
        f.write(content)
    
    return content
    