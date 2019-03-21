#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 30 20:09:24 2018

@author: jkl
"""

import numpy as np
import os

def write2in(param, dirname, typ):
    """ Write a .in file with the `param` specified in the `dirname` for `typ`.
    
    Parameters
    ----------
    param : dict
        Dictionnary of parameters to be used.
    dirname : str
        Path of the working directory to write the .in file.
    typ : str
        Type of file either `R2`, `cR2`, `R3t`, `cR3t`.
    
    Returns
    -------
    String to be writting ot the file.
    """
    if typ == 'R2' or typ == 'cR2':
        write2Din(param, dirname, typ)
    else:
        write3Din(param, dirname, typ)


def write2Din(param, dirname, typ='R2'):
    """ Write a .in file with the `param` specified in the `dirname` for 2D surveys.
    See `write2in` for detailed documentation.
    """
    
    
    # check
    if 'mesh' not in param:
        raise Exception('Need mesh to write configuration file')
    
    # default
    dparam = {
            'lineTitle':'My beautiful survey',
            'job_type':1,
            'mesh_type':4, # meshx, meshy, topo should be defined
            'flux_type':3,
            'singular_type':0,
            'res_matrix':1,
            'scale':1, # only used for triangular mesh
            'num_regions':1,
            'regions':None, # should be defined by the number of element in the mesh
            'patch_x':1,
            'patch_y':1,
#            'nx1':1,
#            'nz1':1,
#            'resistst':1, # default in R2.createMesh()
#            'phasest':-0.1, # default in R2.createMesh()
            'inverse_type':1,
            'target_decrease':0,    
            'data_type':1,
            'reg_mode':0,
            'tolerance':1,
            'max_iter':10,
            'error_mod':2,
            'alpha_aniso':1,
            'alpha_s':10, # penalizing from starting model
            'min_error':0.01, # for IP only
            'a_wgt':0.01, # 0.02 for IP
            'b_wgt':0.02, # 2 for IP
#            'c_wgt':1,
#            'd_wgt':2,
            'rho_min':-1e10,
            'rho_max':1e10,
            'num_xy_poly':-1,
            'xy_poly_table':np.zeros((5,2)),
            'num_elec':None, #should be define when importing data
            'elec_node':None # should be define when building the mesh
            }
    
    
    # check if values are missing
    for a in dparam:
        if a not in param: # parameter missing
            param[a] = dparam[a]
    
    # create text for R2.in
    content = ''
    content = content + '{}\n\n'.format(param['lineTitle'])
    if typ == 'R2':
        content = content + '{}\t{}\t{}\t{}\t{}\n\n'.format(
                            param['job_type'],
                            param['mesh_type'],
                            param['flux_type'],
                            param['singular_type'],
                            param['res_matrix'])
    elif typ == 'cR2':
        content = content + '{}\t{}\t{}\n\n'.format(
                        param['job_type'],
                        param['mesh_type'],
                        param['flux_type'])
    else:
        print('NOT IMPLEMENTED')
    if param['mesh_type'] == 4: # quadrilateral mesh
        meshx = param['meshx']
        meshy = param['meshy']
        topo = param['topo']
        content = content + '\t{}\t{}\t<< numnp_x, numnp_y\n\n'.format(
                len(meshx), len(meshy))
        content = content + ' '.join(['{:.4f}']*len(meshx)).format(*meshx) + '\t<< xx \n\n'
        content = content + ' '.join(['{:.4f}']*len(topo)).format(*topo) + '\t<< topo \n\n'        
        content = content + ' '.join(['{:.4f}']*len(meshy)).format(*meshy) + '\t<< yy \n\n'
    elif param['mesh_type'] == 3:
#        if typ == 'R2':
        content = content + '{}  << scale for triangular mesh\n\n'.format(param['scale'])
    else:
        print('NOT IMPLEMENTED')
    content = content + '{} << num_regions\n'.format(param['num_regions'])
    if param['num_regions'] == 0:
        content = content + param['res0File'] + '\n\n'
    if param['regions'] is None:
        if (param['mesh_type'] == 4) & (typ == 'R2'):
            param['regions'] = np.array([[1, (len(meshx)-1)*(len(meshy)-1), 50]])
        if (param['mesh_type'] == 4) & (typ == 'cR2'):
            param['regions'] = np.array([[1, len(meshx)-1, 1, len(meshy)-1, 1, -0.1]])
    if param['num_regions'] > 0:
        if typ == 'R2':            
            content = content + ''.join(['\t{:.0f}\t{:.0f}\t{} << elem_1, elem_2, value\n']*param['regions'].shape[0]).format(*param['regions'].flatten())
        elif typ =='cR2':
            if param['mesh_type'] == 4:
                content = content + ''.join(['\t{:.0f}\t{:.0f}\t{:.0f}\t{:.0f}\t{}\t{} << nx1, nx2, nz1, nz2, resis, phase\n']*param['regions'].shape[0]).format(*param['regions'].flatten())
            if param['mesh_type'] == 3:
                content = content + ''.join(['\t{:.0f}\t{:.0f}\t{}\t{} << elem_1, elem_2, resist, phase\n']*param['regions'].shape[0]).format(*param['regions'].flatten())
    if (param['mesh_type'] == 4) & (param['job_type'] == 1):
        content = content + '{}\t{}\t<< no. patches in x, no. patches in z\n'.format(param['patch_x'], param['patch_y'])
    if param['job_type'] == 1:
        if param['mesh_type'] == 4|5:
            content = content + '\t{}\t{}\t<< no. patches in x, no. patches in z\n\n'.format(
                    param['patchx'], param['patchy'])
#        if typ == 'R2':
        content = content + '{}\t{}\t<< inverse_type, target_decrease\n\n'.format(
                param['inverse_type'],
                param['target_decrease'])
        if typ == 'R2':
            content = content + '{}\t{}\t<< data type (0=normal;1=log), regularization type\n\n'.format(
                    param['data_type'],
                    param['reg_mode'])
#        elif typ == 'cR2':
#            content = content + '{}\t<< inverse_type\n\n'.format(param['inverse_type'])
        if param['reg_mode'] == 1:
            content = content + '{}\t{}\t{}\t{}\t{}\t<< tolerance, max_iterations, error_mod, alpha_aniso, alpha_s\n\n'.format(
                    param['tolerance'],
                    param['max_iter'],
                    param['error_mod'],
                    param['alpha_aniso'],
                    param['alpha_s'])
        else:
            content = content + '{}\t{}\t{}\t{}\t<< tolerance, max_iterations, error_mod, alpha_aniso\n\n'.format(
                    param['tolerance'],
                    param['max_iter'],
                    param['error_mod'],
                    param['alpha_aniso'])
        if typ == 'R2':
            content = content + '{}\t{}\t{}\t{}\t<<  a_wgt, b_wgt, rho_min, rho_max\n\n'.format(
                    param['a_wgt'],
                    param['b_wgt'],
                    param['rho_min'],
                    param['rho_max'])
        elif typ == 'cR2':
            content = content + '{}\t{}\t{}\t{}\t{}\t<<  min_error, a_wgt, b_wgt, rho_min, rho_max\n\n'.format(
                    param['min_error'],
                    param['a_wgt'],
                    param['b_wgt'],
                    param['rho_min'],
                    param['rho_max'])
#        elif typ == 'cR2':
#                content = content + '{}\t{}\t{}\t{}\t{}\t{}\t<<  a_wgt, b_wgt, c_wgt, d_wgt, rho_min, rho_max\n\n'.format(
#                param['a_wgt'],
#                param['b_wgt'],
#                param['c_wgt'],
#                param['d_wgt'],
#                param['rho_min'],
#                param['rho_max'])
                
    # define polyline
#    if typ == 'R2':
    if param['num_xy_poly'] == -1:
        leftx = param['meshx'][param['node_elec'][0,1]-4]
        rightx = param['meshx'][param['node_elec'][-1,1]+4]
        lefty = np.max(param['topo'])
        righty = lefty
        dist = np.sqrt((rightx-leftx)**2+(righty-lefty)**2)/2
        bottom = np.min(param['topo']) - dist
        param['xy_poly_table'] = np.array([[leftx, lefty],
                                          [rightx, righty],
                                          [rightx, bottom],
                                          [leftx, bottom],
                                          [leftx, lefty]])
        param['num_xy_poly'] = param['xy_poly_table'].shape[0]
    content = content + '{}\t<< num_poly\n'.format(param['num_xy_poly'])
    if param['num_xy_poly'] != 0:
        content = content + ''.join(['{}\t{}\n']*len(param['xy_poly_table'])).format(
                *param['xy_poly_table'].flatten())
        
    # define nodes for electrodes
    param['num_elec'] = param['node_elec'].shape[0]
    content = content + '\n{}\t<< num_electrodes\n'.format(param['num_elec'])
    if param['mesh_type'] == 4:
        content = content + ''.join(['{}\t{}\t{}\n']*len(param['node_elec'])).format(
                *param['node_elec'].flatten())
    elif param['mesh_type'] == 3:
        content = content + ''.join(['{}\t{}\n']*len(param['node_elec'])).format(
                *param['node_elec'].flatten())
    content = content + '\n'
    
    # write configuration file
    fname = typ + '.in'
    with open(os.path.join(dirname,fname),'w') as f:
        f.write(content)
    
    return content



#%%
def write3Din(param, dirname, typ='R3t'):
    """ Write a .in file with the `param` specified in the `dirname` for 3D survey.
    See `write2in` for detailed documentation.
    """

    # default
    dparam = {
            'lineTitle':'My beautiful 3D survey',
            'job_type':1, # inversion by default
            'singular_type':0,
            'num_regions':1,
            'resis':100.0, # default resistivity for regions Ohm.m
            'phase':-2, #default phase for regions mrad
            'inverse_type':0, # different for 3D (it's more the regularization mode)
            'data_type':1,
            'tolerance':1.0,
            'no_improve':1.0 , # 3D specific
            'max_iter':10,
            'error_mod':2,
            'alpha_aniso':1,
            'alpha_s':1, # for 3D inversion_type=2 (difference inversion)
            'cginv_tolerance':0.0001, # 3D specific
            'cginv_maxits':500, # 3D specific
            'alpha_max':1.0e10, # 3D specific
            'num_alpha_steps':10, # 3D specific
            'min_step':0.001, # 3D specific
            'a_wgt':0.01,
            'b_wgt':0.02,
            'c_wgt':1, # cR3t specific
            'rho_min': -1e10, # cR3t specific
            'rho_max': 1e10, # cR3t specific
            'zmin':-10, # 3D specific
            'zmax':0, # 3D specific
            'num_xy_poly':0,
            'elec_node':None # should be define when building the mesh
            }
    
    # check if values are missing
    for a in dparam:
        if a not in param: # parameter missing
            param[a] = dparam[a]
                
    content = ''
    content = content + '{}\n\n'.format(param['lineTitle'])
    content = content + '{}\t{}\t<< job_type, singularity_type\n\n'.format(
            param['job_type'], param['singular_type'])

    # parameters specific to inversion
    if typ=='R3t':
        content = content + '{}\t<< num_regions\n'.format(param['num_regions'])
        
        if param['num_regions'] == 0:
            content = content + param['res0File'] + '\n\n'
        else:
            content = content + '{:f}\n\n'.format(param['resis']) # replace by default rho

        content = content + '{}\t << data_type\n'.format(param['data_type'])
    
        if param['job_type'] == 1:
            content = content + '{}\t<< inverse_type\n\n'.format(param['inverse_type'])
        
            # type of inversion
            if param['inverse_type'] == 0: # normal regularization
                content = content + '{}\t{}\t{}\t{}\t{}\t<< tolerance, no_improve, \
                    max_iterations, error_mod, alpha_aniso\n\n'.format(
                    param['tolerance'], param['no_improve'], param['max_iter'],
                    param['error_mod'], param['alpha_aniso'])
            else: # background regularization or difference inversion
                content = content + '{}\t{}\t{}\t{}\t{}\t{}\t<< tolerance, no_improve, \
                max_iterations, error_mod, alpha_aniso, alpha_s\n\n'.format(
                param['tolerance'], param['no_improve'], param['max_iter'],
                param['error_mod'], param['alpha_aniso'], param['alpha_s'])
                
            content = content + '{}\t{}\t<< cginv_tolerance, cginv_maxits\n\n'.format(
                    param['cginv_tolerance'], param['cginv_maxits'])
            content = content + '{}\t{}\t<< alpha_max, num_alpha_steps\n\n'.format(
                    param['alpha_max'], param['num_alpha_steps'])
            content = content + '{}\t<< min_step\n\n'.format(param['min_step'])
            content = content + '{}\t{}\t<<a_wgt, b_wgt\n\n'.format(param['a_wgt'],
                                 param['b_wgt'])
            content = content + '{}\t{}\t<<zmin, zmax\n'.format(param['zmin'],
                                 param['zmax'])
            
                        
            # define polyline
            content = content + '{}\t<< num_poly\n'.format(param['num_xy_poly'])
            if param['num_xy_poly'] != 0:
                content = content + ''.join(['{}\t{}\n']*len(param['xy_poly_table'])).format(
                        *param['xy_poly_table'].flatten())

    elif typ == 'cR3t':
        content = content + '{}\t<< num_regions\n'.format(param['num_regions'])
        
        if param['num_regions'] == 0:
            content = content + param['res0File'] + '\n\n'
        else:
            content = content + '{}\t{}\t<< resis, phase\n\n'.format(param['resis'], param['phase'] )  # replace by default rho
            
        if param['job_type'] == 1:
            # note that there is no data type in inverse, all is transformed to log by default
            content = content + '{}\t<< inverse_type\n\n'.format(param['inverse_type'])
    
            # type of inversion
            if param['inverse_type'] == 0: # normal regularization
                content = content + '{}\t{}\t{}\t{}\t{}\t<< tolerance, no_improve, \
                    max_iterations, error_mod, alpha_aniso\n\n'.format(
                    param['tolerance'], param['no_improve'], param['max_iter'],
                    param['error_mod'], param['alpha_aniso'])
            else: # background regularization or difference inversion
                content = content + '{}\t{}\t{}\t{}\t{}\t{}\t<< tolerance, no_improve, \
                max_iterations, error_mod, alpha_aniso, alpha_s\n\n'.format(
                param['tolerance'], param['no_improve'], param['max_iter'],
                param['error_mod'], param['alpha_aniso'], param['alpha_s'])
            
            content = content + '{}\t{}\t<< cginv_tolerance, cginv_maxits\n\n'.format(
                    param['cginv_tolerance'], param['cginv_maxits'])
            content = content + '{}\t{}\t<< alpha_max, num_alpha_steps\n\n'.format(
                    param['alpha_max'], param['num_alpha_steps'])
            content = content + '{}\t<< min_step\n\n'.format(param['min_step'])
            
            # cR3t specific
            content = content + '{}\t{}\t{}\t{}\t{}\t<<a_wgt, b_wgt, c_wgt, rho_min, rho_max\n\n'.format(
                            param['a_wgt'], param['b_wgt'], param['c_wgt'],
                            param['rho_min'], param['rho_max'])
            
            content = content + '{}\t{}\t<<zmin, zmax\n'.format(param['zmin'],
                                 param['zmax'])
            
            # define polyline
            content = content + '{}\t<< num_poly\n'.format(param['num_xy_poly'])
            if param['num_xy_poly'] != 0:
                content = content + ''.join(['{}\t{}\n']*len(param['xy_poly_table'])).format(
                        *param['xy_poly_table'].flatten())
    
    # define nodes for electrodes
    param['num_elec'] = param['node_elec'].shape[0]
    content = content + '\n{}\t<< num_electrodes\n'.format(param['num_elec'])
    content = content + ''.join(['1\t{}\t{}\n']*len(param['node_elec'])).format(
            *param['node_elec'].flatten())
    content = content + '\n'
    # write configuration file
    fname = typ + '.in'
    with open(os.path.join(dirname,fname),'w') as f:
        f.write(content)
            
    return content

#test code
#param = {}
#param['node_elec'] = np.c_[np.arange(24)+1, np.arange(24)+1] # needed
#param['inversion_type'] = 0 # test background regularization
#param['job_type'] = 1 # test forward model
#param['num_regions'] = 0
#param['res0File'] = 'resistivity.dat'
#content = write3Din(param,'.')
#content = write3Din(param,'.', typ = 'cR3t') #test cR3.in creation


#%% forward modelling config file 
#def write2in_forward(param, dirname, typ='R2'):
#    """
#    #write a configuration file for a forward model (future update might put this all into one function??)
#    Parameters
#    ----------
#    param : dictionary
#    keys refer to specific settings in R2.exe. The user must provide the number of electrodes
#    and their corresponding nodes. 
#    
#    dirname: string
#    working directory to which the configuration file is saved to
#    
#    typ: string
#    references the which .exe to run. R2 (for resistivity inversion) or cR2 (for
#    induced polarisation inversion)
#    
#    Returns
#    ----------
#    content : string
#    raw text of the configuration file
#    
#    R2.in/cR2.in : file written in directory as flagged by dirname. 
#    
#    """
#    if typ == 'cR2':
#        print('Currently only R2 is supported for write2in_forward')
#    
#    dparam = {
#            'lineTitle':'My beautiful synthetic survey',
#            'job_type':0,
#            'file_name':'_res.dat',#this is the file with the forward model resistivities 
#            'mesh_type':3, # meshx, meshy, topo should be defined
#            'flux_type':3,
#            'singular_type':0,
#            'res_matrix':0,
#            'scale':1, # only used for triangular mesh
#            'num_regions':1,
#            'regions':0, # should be defined by the number of element in the mesh
#            'patch_x':1,
#            'patch_y':1, 
#            'data_type':0,
#            'reg_mode':0,
#            'error_mod':2,
#            'alpha_aniso':1,
#            'num_xy_poly':-1,
#            'xy_poly_table':np.zeros((5,2)),
#            'elec_num':None, #should be define when importing data
#            'elec_node':None, # should be define when building the mesh
#            }
#    
#    for a in dparam:
#        if a not in param: # parameter missing
#            param[a] = dparam[a]
#            
#    if param['job_type']==1:
#        raise Exception ("Use write2in for inverse job r2 configuration files")
#        
#    content = ''
#    content += '{}\n\n'.format(param['lineTitle'])
#    
#    content += '\t{}\t{}\t{:2.1f}\t{}\t{}\n\n'.format(
#                            param['job_type'],
#                            param['mesh_type'],
#                            param['flux_type'],
#                            param['singular_type'],
#                            param['res_matrix'])
#    
#    if param['mesh_type'] == 3:
#        content += '\t{} << vessel_dia\n\n'.format(param['scale'])
#    
#    #input resistivity values to the mesh using the regions tag     
#    if param['num_regions'] >= 1:
#        if typ == 'R2': #code stolen from write2in           
#            content = content + ''.join(['\t{}\t{}\t{} << elem_1, elem_2, value\n']*param['regions'].shape[0]).format(*param['regions'].flatten())
#    elif param['num_regions']==0:
#        content += '\t{} << num_regions\n'.format(param['num_regions'])
#        content += '{:15s}\n\n'.format(param['regionFile'])
#    
#    #insert section here to allow for truncation of the returned mesh? 
#    content += '    0  << num_poly\n\n'
#    
#    #now add the electrodes to the file 
#    content += ' {}  << num_electrodes\n'.format(param['elec_num'])
#    
#    for i in range(len(param['elec_node'])):
#        content += '{}\t{}\n'.format(i+1,param['elec_node'][i])
#    
#    # write configuration file
#    if typ == 'R2':
#        fname = 'R2.in'
#    elif typ == 'cR2':
#        fname = 'cR2.in'
#    with open(os.path.join(dirname,fname),'w') as f:
#        f.write(content)
#        f.close()
#    
#    return content

#%% test code
#import meshTools as mt
#elec_x = np.arange(10)
#elec_y = np.zeros(len(elec_x))
#mesh,meshx,meshy,topo,e_nodes = mt.quad_mesh(elec_x,elec_y)
#param = {}
#param['meshx'] = meshx
#param['meshy'] = meshy
#param['topo'] = topo
#param['mesh_type'] = 4
#param['mesh'] = mesh
#param['node_elec'] = np.c_[1+np.arange(len(e_nodes)), e_nodes, np.ones(len(e_nodes))].astype(int)
#content = write2in(param=param, dirname='.')


