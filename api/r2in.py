#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 30 20:09:24 2018

@author: jkl
"""

import numpy as np
import os

def write2in(param, dirname, typ='R2'):
    ''' write R2.in file
    param = dictionnary of the parameters
    '''
    
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
            'inverse_type':1,
            'target_decrease':1,    
            'data_type':0,
            'reg_mode':0,
            'tol':1,
            'max_iter':10,
            'error_mod':2,
            'alpha_aniso':1,
            'a_wgt':0.01,
            'b_wgt':0.02,
            'd_wgt':2,
            'c_wgt':1,
            'rho_min':-1000,
            'rho_max':1000,
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
        content = content + '{}  << scale for triangular mesh\n\n'.format(param['scale'])
    else:
        print('NOT IMPLEMENTED')
    content = content + '{} << num_regions\n'.format(param['num_regions'])
    if param['regions'] is None:
        if param['mesh_type'] == 4:
            param['regions'] = np.array([[1, (len(meshx)-1)*(len(meshy)-1), 50]])
    if param['num_regions'] > 0:
        content = content + ''.join(['\t{}\t{}\t{} << elem_1, elem_2, value\n']*param['regions'].shape[0]).format(*param['regions'].flatten())
    content = content + '{}\t{}\t<< no. patches in x, no. patches in z\n'.format(param['patch_x'], param['patch_y'])
    if param['job_type'] == 1 & param['mesh_type'] == 4|5:
        content = content + '\t{}\t{}\t<< no. patches in x, no. patches in z\n\n'.format(
                param['patchx'], param['patchy'])
    
    content = content + '{}\t{}\t<< inverse_type, target_decrease\n\n'.format(
            param['inverse_type'],
            param['target_decrease'])    
    content = content + '{}\t{}\t<< data type (0=normal;1=log), regularization type\n\n'.format(
            param['data_type'],
            param['reg_mode'])
    content = content + '{}\t{}\t{}\t{}\t<< tolerance, max_iterations, error_mod, alpha_aniso\n\n'.format(
            param['tol'],
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
            content = content + '{}\t{}\t{}\t{}\t{}\t{}\t<<  a_wgt, b_wgt, c_wgt, d_wgt, rho_min, rho_max\n\n'.format(
            param['a_wgt'],
            param['b_wgt'],
            param['c_wgt'],
            param['d_wgt'],
            param['rho_min'],
            param['rho_max'])
    if typ == 'R2':
        if param['num_xy_poly'] == -1:
            leftx = param['meshx'][param['node_elec'][0,1]-4]
            rightx = param['meshx'][param['node_elec'][-1,1]+4]
            lefty = param['topo'][param['node_elec'][0,1]-4]
            righty = param['topo'][param['node_elec'][-1,1]+4]
            dist = np.sqrt((rightx-leftx)**2+(righty-lefty)**2)/2
            param['xy_poly_table'] = np.array([[leftx, lefty],
                                              [rightx, righty],
                                              [rightx, -dist],
                                              [leftx, -dist],
                                              [leftx, lefty]])
            param['num_xy_poly'] = param['xy_poly_table'].shape[0]
        content = content + '{}\t<< num_poly\n'.format(param['num_xy_poly'])
        if param['num_xy_poly'] != 0:
            content = content + ''.join(['{}\t{}\n']*len(param['xy_poly_table'])).format(
                    *param['xy_poly_table'].flatten())
    param['num_elec'] = param['node_elec'].shape[0]
    content = content + '\n{}\t<< num_electrodes\n'.format(param['num_elec'])
    content = content + ''.join(['{}\t{}\t{}\n']*len(param['node_elec'])).format(
            *param['node_elec'].flatten())
    content = content + '\n'
    

    # write configuration file
    if typ == 'R2':
        fname = 'R2.in'
    elif typ == 'cR2':
        fname = 'cR2.in'
    with open(os.path.join(dirname,fname),'w') as f:
        f.write(content)
    
    return content


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


