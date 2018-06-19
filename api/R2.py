#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 30 16:48:54 2018

@author: jkl
"""

import numpy as np
import pandas as pd
import os
import sys
import shutil
import platform
import matplotlib.pyplot as plt
from subprocess import PIPE, call, Popen
OS = platform.system()

#sys.path.append(os.path.relpath('../api'))

#from QuadMesh import QuadMesh
#from TriMesh import TriMesh

from Survey import Survey
from r2in import write2in
import meshTools as mt
from meshTools import Mesh_obj
from gmshWrap import tri_mesh
#from testMesh import QuadMesh



class R2(object): # R2 master class instanciated by the GUI
    def __init__(self, dirname=''):
        if dirname == '':
            dirname = os.getcwd()
            print('using the current directory:', dirname)
        self.dirname = dirname # working directory
        self.surveys = [] # list of survey object
        self.surveysInfo = [] # info about surveys (date)
        self.mesh = None # mesh object (one per R2 instance)
        self.param = {} # dict configuration variables for inversion
        self.configFile = ''
        self.typ = 'R2' # or cR2 or R3, cR3
        self.errTyp = 'obs' # type of error to add for DC
        self.errTypIP = 'none' # type of error to add for IP phase
        
        
    def setwd(self, dirname):
        ''' set the working directory
        '''
        self.dirname = dirname
    
    
    def createSurvey(self, fname, ftype='Syscal', info={}):
        ''' read electrode and quadrupoles data and return 
        a survey object
        
        fname : filename to be parsed
        ftype : type of file to be parsed
        info : dict of info about the survey
        '''    
        self.surveys.append(Survey(fname, ftype))
        self.surveysInfo.append(info)
        
        # define electrode position according to first survey
        if len(self.surveys) == 1:
            self.elec = self.surveys[0].elec
            
            # attribute method of Survey object to R2
            self.pseudo = self.surveys[0].pseudo
            self.plotError = self.surveys[0].plotError
            self.linfit = self.surveys[0].linfit
            self.lmefit = self.surveys[0].lmefit
            self.pwlfit = self.surveys[0].pwlfit
            self.phaseplotError = self.surveys[0].phaseplotError
            self.plotIPFit = self.surveys[0].plotIPFit
        
    
    def createMesh(self, typ='default', **kwargs):
        ''' create a mesh object
        typ:
            quad : quadrilateral mesh
            triang : triangular mesh
        '''
        if typ == 'default':
            if self.elec[:,2].sum() == 0:
                typ = 'quad'
                print('Using a quadrilateral mesh')
            else:
                typ = 'trian'
                print('Using a triangular mesh')
        if typ == 'quad':
#            mesh = QuadMesh(elec, nnode=4)
            elec_x = self.elec[:,0]
            elec_y = self.elec[:,1]
            mesh,meshx,meshy,topo,e_nodes = mt.quad_mesh(elec_x,elec_y,elemx=8)
#            mesh = QuadMesh()
#            meshx, meshy, topo, e_nodes = mesh.createMesh(elec=self.elec, **kwargs)            
            self.param['meshx'] = meshx
            self.param['meshy'] = meshy
            self.param['topo'] = topo
            self.param['mesh_type'] = 4
            self.param['node_elec'] = np.c_[1+np.arange(len(e_nodes)), e_nodes, np.ones(len(e_nodes))].astype(int)
        if typ == 'trian':
            mesh, e_ranges = tri_mesh([],[], self.elec[:,0], self.elec[:,1], save_path=os.path.join(self.dirname, 'mesh.dat'))
            self.param['mesh_type'] = 3
            self.param['num_regions'] = len(e_ranges)
            self.param['regions'] = np.array(e_ranges)
            self.param['num_xy_poly'] = 0
            e_nodes = np.arange(len(self.elec))+1
            self.param['node_elec'] = np.c_[1+np.arange(len(e_nodes)), e_nodes, np.ones(len(e_nodes))].astype(int)
        self.mesh = mesh
        self.param['mesh'] = mesh
        
        
    def showMesh(self, ax=None):
        if self.mesh is None:
            raise Exception('Mesh undefined')
        else:
#            xlim = (np.min(self.elec[:,0]-20, np.max(self.elec[:,0])))
#            ylim = (0, 110) # TODO
#            self.mesh.show(xlim=xlim, ylim=ylim) # add ax argument
            self.mesh.show(ax=ax)
    
    def write2in(self, param={}, typ=''):
        ''' create configuration file for inversion
        '''
        if typ == '':
            typ = self.typ
        # check if survey has reciprocal (and so error model)
        if all(self.surveys[0].df['irecip'].values == 0):
            self.param['wgt_a'] = 0.01
            self.param['wgt_b'] = 0.02
        if typ == 'cR2':
            if self.errTypIP != 'none': # we have individual errors
                self.param['c_wgt'] = 0
                self.param['d_wgt'] = 0
                self.param['a_wgt'] = 0.01 # not sure of that (Gui)
            else:
                print('kkk')
                self.param['c_wgt'] = 0.1 # better if set by user !!
                self.param['d_wgt'] = 0.2
        # all those parameters are default but the user can change them and call
        # write2in again
        for p in param:
            self.param[p] = param[p]
            
        self.configFile = write2in(self.param, self.dirname, typ=typ)


    def write2protocol(self, **kwargs): # DEPRACATED
        self.surveys[0].write2protocol(**kwargs)
        
        
    def runR2(self):
        # run R2.exe
        exeName = self.typ + '.exe'
        cwd = os.getcwd()
        os.chdir(self.dirname)
        targetName = os.path.join(self.dirname, exeName)
        
        # copy R2.exe
        if ~os.path.exists(targetName):
            shutil.copy(os.path.join(cwd, 'exe', exeName), targetName)  
        
        if OS == 'Linux':
            cmd = ['wine',exeName]
        if OS == 'Windows':
            cmd = [exeName]
        
        p = Popen(cmd, stdout=PIPE, shell=False)
        while p.poll() is None:
            line = p.stdout.readline().rstrip()
            print(line.decode('utf-8'))
            
        os.chdir(cwd)
        
        
    def invert(self, param={}, iplot=True):
        ''' invert the data, first generate R2.in file, then run
        inversion using appropriate wrapper, then return results
        '''
        # create mesh if not already done
        if 'mesh' not in self.param:
            self.createMesh()
        
        # write configuration file
        if self.configFile == '':
            self.write2in(param=param)
        
        if self.typ == 'R2':
            ipBool = False
        elif self.typ == 'cR2':
            ipBool = True
        else:
            print('NOT IMPLEMENTED YET')
            
        self.surveys[0].write2protocol(os.path.join(self.dirname, 'protocol.dat'),
                errTyp=self.errTyp, ip=ipBool, errTypIP=self.errTypIP)
            
              
        self.runR2()
        
        if iplot is True:
#            self.showResults()
            self.showSection()
    
    
    def showResults(self, ax=None, **kwargs):
        fresults = os.path.join(self.dirname, 'f001_res.vtk')
        if os.path.isfile(fresults):
            mesh_dict=mt.vtk_import(fresults)#makes a dictionary of a mesh 
            mesh = Mesh_obj.mesh_dict2obj(mesh_dict)# this is a mesh_obj class instance 
            mesh.show(ax=ax, **kwargs)
        else:
            print('Sorry no VTK output produced')
        
    def showSection(self, fname='', ax=None, ilog10=True, isen=False, figsize=(8,3)):
        if fname == '':
            fname = os.path.join(self.dirname, 'f001_res.dat')
        res = pd.read_csv(fname, delimiter=' *', header=None, engine='python').values
        lenx = len(np.unique(res[:,0]))
        leny = len(np.unique(res[:,1]))
        x = res[:,0].reshape((leny, lenx), order='F')
        y = res[:,1].reshape((leny, lenx), order='F')
        z = res[:,2].reshape((leny, lenx), order='F')
        if isen:
            sen = pd.read_csv(fname.replace('res','sen'), delimiter=' *', header=None, engine='python').values
            lenx = len(np.unique(sen[:,0]))
            leny = len(np.unique(sen[:,1]))
        #            xs = sen[:,0].reshape((leny, lenx), order='F')
        #            ys = sen[:,1].reshape((leny, lenx), order='F')
            zs = sen[:,2].reshape((leny, lenx), order='F')
            zs = np.log10(zs)
            zs -= np.min(zs)
            alpha = zs/np.max(zs)
        #            alpha[alpha < 0] = 0
            print(np.max(alpha), np.min(alpha))
        if ilog10:
            z = np.log10(z)
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.get_figure()
        cax = ax.pcolormesh(x, y, z)
        ax.plot(self.elec[:,0], self.elec[:,1], 'ko')
    #    fig.canvas.draw() # need to draw the figure to have the cax.get_facecolors()
    #    print(cax.get_facecolors().shape)
    #    print(alpha.flatten().shape)
    #    for a in cax.get_facecolors():
    #        a[3] = 0
        #for a, b in zip(cax.get_facecolors(), alpha.flatten()):
        #    a[3] = 0.5
        #    print(a)
    #    fig.canvas.draw()
        cbar = fig.colorbar(cax, ax=ax)
        if ilog10:
            cbar.set_label(r'$\log_{10}(\rho) [\Omega.m]$')
        else:
            cbar.set_label(r'$\rho [\Omega.m]$')
        ax.set_ylabel('Depth [m]')
        ax.set_xlabel('Distance [m]')
#        fig.tight_layout()
    #    fig.show()
#        return fig
    
    def pseudoError(self, ax=None):
        ''' plot pseudo section of errors from file f001_err.dat
        '''
        err = np.genfromtxt(os.path.join(self.dirname, 'f001_err.dat'), skip_header=1)
        array = err[:,[-2,-1,-4,-3]].astype(int)
        errors = err[:,0]
        spacing = np.diff(self.elec[[0,1],0])
        pseudo(array, errors, spacing, ax=ax, label='Normalized Errors', log=False, geom=False, contour=False)


def pseudo(array, resist, spacing, label='', ax=None, contour=False, log=True, geom=True):
    nelec = np.max(array)
    elecpos = np.arange(0, spacing*nelec, spacing)
    resist = resist
    
    if geom: # compute and applied geometric factor
        apos = elecpos[array[:,0]-1]
        bpos = elecpos[array[:,1]-1]
        mpos = elecpos[array[:,2]-1]
        npos = elecpos[array[:,3]-1]
        AM = np.abs(apos-mpos)
        BM = np.abs(bpos-mpos)
        AN = np.abs(apos-npos)
        BN = np.abs(bpos-npos)
        K = 2*np.pi/((1/AM)-(1/BM)-(1/AN)+(1/BN)) # geometric factor
        resist = resist*K
        
    if log:
        resist = np.sign(resist)*np.log10(np.abs(resist))
    if label == '':
        if log:
            label = r'$\log_{10}(\rho_a)$ [$\Omega.m$]'
        else:
            label = r'$\rho_a$ [$\Omega.m$]'
    
    cmiddle = np.min([elecpos[array[:,0]-1], elecpos[array[:,1]-1]], axis=0) \
        + np.abs(elecpos[array[:,0]-1]-elecpos[array[:,1]-1])/2
    pmiddle = np.min([elecpos[array[:,2]-1], elecpos[array[:,3]-1]], axis=0) \
        + np.abs(elecpos[array[:,2]-1]-elecpos[array[:,3]-1])/2
    xpos = np.min([cmiddle, pmiddle], axis=0) + np.abs(cmiddle-pmiddle)/2
    ypos = - np.sqrt(2)/2*np.abs(cmiddle-pmiddle)
    
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure
    cax = ax.scatter(xpos, ypos, c=resist, s=70)#, norm=mpl.colors.LogNorm())
    cbar = fig.colorbar(cax, ax=ax)
    cbar.set_label(label)
    ax.set_title('Pseudo Section')

        
#%% test code
#k = R2('/media/jkl/data/phd/tmp/r2gui/api/test')
#k.typ = 'cR2'
#k.createSurvey('test/syscalFile.csv', ftype='Syscal')
#k.createSurvey('test/rifleday8_n2.csv', ftype='Syscal')
#k.pseudo(contour=True)
#k.linfit(iplot=True)
#k.lmefit(iplot=True)
#k.createMesh(typ='quad')
#k.createMesh(typ='trian')
#k.mesh.show()
#fig, ax = plt.subplots()
#fig.suptitle('kkk')
#k.mesh.show(ax=ax)
#k.write2in()
#k.invert(iplot=False)
#k.pseudoError()
#k.showSection()
#fig, ax = plt.subplots()
#fig.suptitle('kkk')
#k.showResults(ax=ax)

