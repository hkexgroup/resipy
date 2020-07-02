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
    
    # default
    dparam = {
            'lineTitle':'My beautiful survey',
            'job_type':1,
            'mesh_type':6, # meshx, meshy, topo should be defined
            'flux_type':3,
            'singular_type':0,
            'res_matrix':1,
            'scale':1, # only used for triangular mesh
            'num_regions':1,
            'regions':None, # should be defined by the number of element in the mesh
            'patch_x':1,
            'patch_z':1,
            'inverse_type':1,
            'target_decrease':0,
            'qual_ratio':0,
            'data_type':1,
            'reg_mode':0,
            'tolerance':1,
            'max_iter':10,
            'error_mod':2,
            'alpha_aniso':1,
            'alpha_s':10, # penalizing from starting model
            'min_error':0.01, # for IP only
            # 'a_wgt':0.01, # 0.02 for IP
            # 'b_wgt':0.02, # 2 for IP
            'rho_min':-1e10,
            'rho_max':1e10,
            'num_xz_poly':-1,
            'xz_poly_table':np.zeros((5,2)),
            'num_elec':None, #should be define when importing data
            'elec_node':None # should be define when building the mesh
            }
    
    
    # check if values are missing
    for a in dparam:
        if a not in param: # parameter missing
            param[a] = dparam[a]
    
    # select default a_wgt/b_wgt (note that the DC a_wgt (offset) does not have
    # equivalent for IP and the the DC b_wgt is the IP a_wgt)
    if 'a_wgt' not in param:
        param['a_wgt'] = 0.01 if typ == 'R2' else 0.02
    if 'b_wgt' not in param:
        param['b_wgt'] = 0.02 if typ == 'R2' else 2
        
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
        meshy = param['meshz']
        topo = param['topo']
        content = content + '\t{}\t{}\t<< numnp_x, numnp_z\n\n'.format(
                len(meshx), len(meshy))
        content = content + ' '.join(['{:.4f}']*len(meshx)).format(*meshx) + '\t<< xx \n\n'
        content = content + ' '.join(['{:.4f}']*len(topo)).format(*topo) + '\t<< topo \n\n'        
        content = content + ' '.join(['{:.4f}']*len(meshy)).format(*meshy) + '\t<< yy \n\n'
    elif param['mesh_type'] == 3: # triangular mesh
        content = content + '{}  << scale for triangular mesh\n\n'.format(param['scale'])
    elif param['mesh_type'] == 6: # general quadrilateral mesh (with mesh.dat)
        pass # no scale for this
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
        content = content + '{}\t{}\t<< inverse_type, target_decrease\n\n'.format(
                param['inverse_type'],
                param['target_decrease'])
        if param['inverse_type'] == 3:
            content = content + '{}\t<< qualitative ratio\n'.format(param['qual_ratio'])
            content = content + '{}\t{}\t<< rho_min, rho_max\n'.format(param['rho_max'], param['rho_min'])
        else:
            if typ == 'R2':
                content = content + '{}\t{}\t<< data type (0=normal;1=log), regularization type\n\n'.format(
                        param['data_type'],
                        param['reg_mode'])
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

    if param['num_xz_poly'] == -1:
        leftx = param['meshx'][param['node_elec'][0,1]-4]
        rightx = param['meshx'][param['node_elec'][-1,1]+4]
        lefty = np.max(param['topo'])
        righty = lefty
        dist = np.sqrt((rightx-leftx)**2+(righty-lefty)**2)/2
        bottom = np.min(param['topo']) - dist
        param['xz_poly_table'] = np.array([[leftx, lefty],
                                          [rightx, righty],
                                          [rightx, bottom],
                                          [leftx, bottom],
                                          [leftx, lefty]])
        param['num_xz_poly'] = param['xz_poly_table'].shape[0]
    content = content + '{}\t<< num_poly\n'.format(param['num_xz_poly'])
    if param['num_xz_poly'] != 0:
        content = content + ''.join(['{}\t{}\n']*len(param['xz_poly_table'])).format(
                *param['xz_poly_table'].flatten())
        
    # define nodes for electrodes
    param['num_elec'] = len(param['node_elec'][0])
    content = content + '\n{}\t<< num_electrodes\n'.format(param['num_elec'])
    # if param['mesh_type'] == 4:
        # content = content + ''.join(['{}\t{}\t{}\n']*len(param['node_elec'])).format(
                # *param['node_elec'].flatten())
    if param['mesh_type'] == 3 or param['mesh_type'] == 6:
        for i in range(param['num_elec']):
            content = content + '{}\t{}\n'.format(param['node_elec'][0][i],
                                                  param['node_elec'][1][i])
    else:
        print('Quadrilateral mesh not supported')
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
            'inverse_type':1, # not used by R3t but for consistency with R2
            'target_decrease':0,
            'data_type':1,
            'reg_mode':0,
            'tolerance':1.0,
            'max_iter':10,
            'error_mod':2,
            'alpha_aniso':1,
            'alpha_s':1, # for 3D inversion_type=2 (difference inversion)
            'min_error':0.01, # for IP only
            # 'a_wgt':0.01, # DC
            # 'b_wgt':0.02, # DC
            # 'a_wgt':0.02, # IP
            # 'b_wgt':2, # IP
            'rho_min': -1e10, # cR3t specific
            'rho_max': 1e10, # cR3t specific
            'zmin':-10, # 3D specific
            'zmax':0, # 3D specific
            'num_xy_poly':0, # 3D specific
            'xy_poly_table':None,
            'elec_node':None # should be define when building the mesh
            }
    
    # check if values are missing
    for a in dparam:
        if a not in param: # parameter missing
            param[a] = dparam[a]
            
    # select default a_wgt/b_wgt (note that the DC a_wgt (offset) does not have
    # equivalent for IP and the the DC b_wgt is the IP a_wgt)
    if 'a_wgt' not in param:
        param['a_wgt'] = 0.01 if typ == 'R2' else 0.02
    if 'b_wgt' not in param:
        param['b_wgt'] = 0.02 if typ == 'R2' else 2
            
                
    content = ''
    content = content + '{}\n\n'.format(param['lineTitle'])
    content = content + '{}\t{}\t<< job_type, singularity_type\n\n'.format(
            param['job_type'], param['singular_type'])

    # parameters specific to inversion
    content = content + '{}\t<< num_regions\n'.format(param['num_regions'])
    
    if param['num_regions'] == 0:
        content = content + param['res0File'] + '\n\n'
    else:
        content = content + '{:f}\n\n'.format(param['resis']) # replace by default rho

    if param['job_type'] == 1:
        content = content + '{}\t{}\t<< inverse_type, target_decrease\n'.format(
            param['inverse_type'], param['target_decrease'])
        if typ == 'R3t':
            content = content + '{}\t{}\t<< data_type, reg_mode\n\n'.format(
                param['data_type'], param['reg_mode'])
    
        # type of inversion
        if param['reg_mode'] == 1: # background regularization
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
        
        # error weighting
        if typ == 'R3t':
            content = content + '{}\t{}\t{}\t{}\t<<  a_wgt, b_wgt, rho_min, rho_max\n\n'.format(
                    param['a_wgt'],
                    param['b_wgt'],
                    param['rho_min'],
                    param['rho_max'])
        elif typ == 'cR3t':
            content = content + '{}\t{}\t{}\t{}\t{}\t<<  min_error, a_wgt, b_wgt, rho_min, rho_max\n\n'.format(
                    param['min_error'],
                    param['a_wgt'],
                    param['b_wgt'],
                    param['rho_min'],
                    param['rho_max'])

        content = content + '{}\t{}\t<<zmin, zmax\n'.format(param['zmin'],
                             param['zmax'])
        
        # define polyline
        content = content + '{}\t<< num_poly\n'.format(param['num_xy_poly'])
        if param['num_xy_poly'] != 0:
            content = content + ''.join(['{}\t{}\n']*len(param['xy_poly_table'])).format(
                    *param['xy_poly_table'].flatten())

    # define nodes for electrodes
    param['num_elec'] = len(param['node_elec'][0])
    content = content + '\n{}\t<< num_electrodes\n'.format(param['num_elec'])
    if len(param['node_elec'][0][0].split(' ')) == 1: # the string number is not part of the label so let's add it
        param['node_elec'][0] = ['1 ' + a for a in param['node_elec'][0]]
    for i in range(param['num_elec']):
        content = content + '{}\t{}\n'.format(param['node_elec'][0][i],
                                              param['node_elec'][1][i])
    content = content + '\n'
    # write configuration file
    fname = typ + '.in'
    with open(os.path.join(dirname,fname),'w') as f:
        f.write(content)
            
    return content


