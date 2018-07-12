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

sys.path.append(os.path.relpath('..'))

from api.Survey import Survey
from api.r2in import write2in
import api.meshTools as mt
from api.meshTools import Mesh_obj
from api.gmshWrap import tri_mesh


class R2(object): # R2 master class instanciated by the GUI
    def __init__(self, dirname=''):
        if dirname == '':
            dirname = os.getcwd()
            print('using the current directory:', dirname)
        self.dirname = dirname # working directory (for the datas)
        self.cwd = os.getcwd() # directory of the code
        self.surveys = [] # list of survey object
        self.surveysInfo = [] # info about surveys (date)
        self.mesh = None # mesh object (one per R2 instance)
        self.param = {} # dict configuration variables for inversion
        self.configFile = ''
        self.typ = 'R2' # or cR2 or R3, cR3
        self.errTyp = '' # type of error to add for DC
        self.errTypIP = 'none' # type of error to add for IP phase
        self.iBorehole = False
        self.iTimeLapse = False
        self.meshResults = [] # contains vtk mesh object of inverted section
        
        
    def setwd(self, dirname):
        ''' set the working directory
        '''
        self.dirname = dirname
    
    
    def createSurvey(self, fname, ftype='Syscal', info={}, spacing=None):
        ''' read electrode and quadrupoles data and return 
        a survey object
        
        fname : filename to be parsed
        ftype : type of file to be parsed
        info : dict of info about the survey
        '''    
        self.surveys.append(Survey(fname, ftype, spacing=spacing))
        self.surveysInfo.append(info)
        
        # define electrode position according to first survey
        if len(self.surveys) == 1:
            self.elec = self.surveys[0].elec
            
            # attribute method of Survey object to R2
            self.pseudoIP = self.surveys[0].pseudoIP
            self.pseudo = self.surveys[0].pseudo
            self.plotError = self.surveys[0].plotError
            self.linfit = self.surveys[0].linfit
            self.lmefit = self.surveys[0].lmefit
            self.pwlfit = self.surveys[0].pwlfit
            self.phaseplotError = self.surveys[0].phaseplotError
            self.plotIPFit = self.surveys[0].plotIPFit
            self.heatmap = self.surveys[0].heatmap

    def createTimeLapseSurvey(self, dirname, ftype='Syscal', info={}, spacing=None, isurveys=[]):
        ''' read electrode and quadrupoles data and return 
        a survey object
        
        fname : filename to be parsed
        ftype : type of file to be parsed
        info : dict of info about the survey
        '''
        self.iTimeLapse = True
        self.iTimeLapseReciprocal = [] # true if survey has reciprocal
        files = np.sort(os.listdir(dirname))
        for f in files:
            self.createSurvey(os.path.join(dirname, f))
            haveReciprocal = all(self.surveys[-1].df['irecip'].values == 0)
            self.iTimeLapseReciprocal.append(haveReciprocal)
            print('---------', f, 'imported')
            if len(self.surveys) == 1:
                ltime = len(self.surveys[0].df)
            if len(self.surveys) > 1:
                if len(self.surveys[-1].df) != ltime:
                    print('ERROR:', f, 'survey doesn\'t have the same length')
                    return
        self.iTimeLapseReciprocal = np.array(self.iTimeLapseReciprocal)
        self.elec = self.surveys[0].elec
        
        
        # create bigSurvey
        print('creating bigSurvey')
        self.bigSurvey = Survey(os.path.join(dirname, files[0]), ftype=ftype)
        # then override the df
        if len(isurveys) == 0: # assume all surveys would be use for error modelling
            isurveys = np.ones(len(self.surveys), dtype=bool)
        isurveys = np.where(isurveys)[0] # convert to indices
        df = self.bigSurvey.df.copy()
        for i in isurveys:
            df.append(self.surveys[i].df)
        self.bigSurvey.df = df.copy() # override it
        self.bigSurvey.dfOrigin = df.copy()
        
        self.pseudo = self.surveys[0].pseudo # just display first pseudo section
            
        self.plotError = self.bigSurvey.plotError
        self.linfit = self.bigSurvey.linfit
        self.lmefit = self.bigSurvey.lmefit
        self.pwlfit = self.bigSurvey.pwlfit
        self.phaseplotError = self.bigSurvey.phaseplotError
        self.plotIPFit = self.bigSurvey.plotIPFit
            
    
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
            mesh,meshx,meshy,topo,e_nodes = mt.quad_mesh(elec_x,elec_y,**kwargs)
#            mesh = QuadMesh()
#            meshx, meshy, topo, e_nodes = mesh.createMesh(elec=self.elec, **kwargs)            
            self.param['meshx'] = meshx
            self.param['meshy'] = meshy
            self.param['topo'] = topo
            self.param['mesh_type'] = 4
            self.param['node_elec'] = np.c_[1+np.arange(len(e_nodes)), e_nodes, np.ones(len(e_nodes))].astype(int)
        if typ == 'trian':
            mesh, e_ranges = tri_mesh([],[], self.elec[:,0], self.elec[:,1], path=os.path.join(self.cwd, 'api', 'exe'), save_path=os.path.join(self.dirname, 'mesh.dat'))
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
            self.mesh.show(ax=ax, color_bar=False)
    
    def write2in(self, param={}, typ=''):
        ''' create configuration file for inversion
        '''
        if typ == '':
            typ = self.typ
        if all(self.surveys[0].df['irecip'].values == 0):
            if 'a_wgt' not in self.param:
                self.param['a_wgt'] = 0.01
            if 'b_wft' not in self.param:
                self.param['b_wgt'] = 0.02
        if typ == 'cR2':
            if self.errTypIP != 'none': # we have individual errors
                if 'b_wgt' not in self.param:
                    self.param['b_wgt'] = 0
                if 'c_wgt' not in self.param:
                    self.param['c_wgt'] = 0
                if 'a_wgt' not in self.param:
                    self.param['a_wgt'] = 0.01 # not sure of that (Gui)
            else:
                if 'c_wgt' not in self.param:
                    self.param['c_wgt'] = 1 # better if set by user !!
                if 'd_wgt' not in self.param:
                    self.param['d_wgt'] = 2
                
        # all those parameters are default but the user can change them and call
        # write2in again
        for p in param:
            self.param[p] = param[p]
        
        if self.iTimeLapse == True:
            refdir = os.path.join(self.dirname, 'ref')
            if os.path.exists(refdir) == False:
                os.mkdir(refdir)
            param = self.param
            param['num_xy_poly'] = 0
            self.configFile = write2in(param, refdir, typ=typ)
            param = self.param
            param['num_regions'] = 0
            param['regularization_type'] = 2
            param['timeLapse'] = 'Start_res.dat'
            write2in(param, self.dirname, typ=typ)
        else:
            self.configFile = write2in(self.param, self.dirname, typ=typ)
        
        

    def write2protocol(self, errTyp='', errTypIP='', **kwargs):
        if self.typ == 'R2':
            ipBool = False
        elif self.typ == 'cR2':
            ipBool = True
        else:
            print('NOT IMPLEMENTED YET')

        if errTyp == '':
            errTyp = self.errTyp
        if ipBool == True:
            if errTypIP == '':
                errTypIP = self.errTypIP
        
        
        if self.iTimeLapse == False:
            self.surveys[0].write2protocol(os.path.join(self.dirname, 'protocol.dat'),
                        errTyp=self.errTyp, ip=ipBool, errTypIP=self.errTypIP)
        else:
            # a bit simplistic but assign error to all based on Transfer resistance
            allHaveReciprocal = all(self.iTimeLapseReciprocal == True)
            # let's assume it's False all the time for now
            content = ''
            for i, s in enumerate(self.surveys):
                content = content + str(len(s.df)) + '\n'
                if errTyp != '': # there is an error model
                    s.df['error'] = self.bigSurvey.errorModel(s.df['resist'].values)
                    s.df['index'] = np.arange(1, len(s.df)+1)
                    content = content + s.df[['index','a','b','m','n','resist','error']].to_csv(sep='\t', header=False, index=False)
                else:
                    s.df['index'] = np.arange(1, len(s.df)+1)
                    content = content + s.df[['index','a','b','m','n','resist']].to_csv(sep='\t', header=False, index=False)
                if i == 0:
                    refdir = os.path.join(self.dirname, 'ref')
                    if os.path.exists(refdir) == False:
                        os.mkdir(refdir)
                    with open(os.path.join(refdir, 'protocol.dat'), 'w') as f:
                        f.write(content) # write the protocol for the reference file
            with open(os.path.join(self.dirname, 'protocol.dat'), 'w') as f:
                f.write(content)
        
        
    def runR2(self, dirname='', dump=print):
        # run R2.exe
        exeName = self.typ + '.exe'
        cwd = os.getcwd()
        if dirname == '':
            dirname = self.dirname
        os.chdir(dirname)
        targetName = os.path.join(dirname, exeName)
        actualPath = os.path.dirname(os.path.relpath(__file__))
        
        # copy R2.exe
        if ~os.path.exists(targetName):
            shutil.copy(os.path.join(actualPath, 'exe', exeName), targetName)  
        
        if OS == 'Linux':
            cmd = ['wine',exeName]
        if OS == 'Windows':
            cmd = [exeName]
        
        p = Popen(cmd, stdout=PIPE, shell=False)
        while p.poll() is None:
            line = p.stdout.readline().rstrip()
            dump(line.decode('utf-8'))
            
        os.chdir(cwd)
        
        
    def invert(self, param={}, iplot=True, dump=print):
        ''' invert the data, first generate R2.in file, then run
        inversion using appropriate wrapper, then return results
        '''
        # create mesh if not already done
        if 'mesh' not in self.param:
            self.createMesh()
        
        # write configuration file
        if self.configFile == '':
            self.write2in(param=param)
        
        self.write2protocol()    
             
        if self.iTimeLapse == True:
            refdir = os.path.join(self.dirname, 'ref')
            self.runR2(refdir, dump=dump)
            print('----------- finished inverting reference model ------------')
            shutil.copy(os.path.join(refdir, 'f001_res.dat'),
                    os.path.join(self.dirname, 'Start_res.dat'))
            self.runR2(dump=dump)
        else:
            self.runR2(dump=dump)
        
        if iplot is True:
#            self.showResults()
            self.showSection() # TODO need to debug that for timelapse and even for normal !
            # pass an index for inverted survey time
    
    def showResults(self, index=0, ax=None, edge_color='none', attr='Resistivity(log10)', sens=True, **kwargs):
        if len(self.meshResults) == 0:
            self.getResults()
        if len(self.meshResults) > 0:
            self.meshResults[index].show(ax=ax, edge_color=edge_color, attr=attr, sens=sens, **kwargs)
        else:
            print('Unexpected Error')

    
    def getResults(self):
        for i in range(100):
            fresults = os.path.join(self.dirname, 'f' + str(i+1).zfill(3) + '_res.vtk')
            if os.path.exists(fresults):
                mesh = mt.vtk_import(fresults)
                mesh.elec_x = self.elec[:,0]
                mesh.elec_y = self.elec[:,1]
                self.meshResults.append(mesh)
            else:
                break

            
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
#k.pwlfit()
#k.errTyp='obs'
#k.lmefit(iplot=True)
#k.createMesh(typ='quad')
#k.createMesh(typ='trian')
#k.mesh.show()
#fig, ax = plt.subplots()
#fig.suptitle('kkk')
#k.mesh.show(ax=ax)
#k.write2in()
#k.plotIPFit()
#k.errTyp = 'pwl'
#k.errTypIP = 'pwl'
#k.invert(iplot=False)
#k.pseudoError()
#k.showSection()
#fig, ax = plt.subplots()
#fig.suptitle('hkk')
#k.showResults()
#k.showResults(edge_color='none', sens=True)
#k.showResults(attr=attr[0])
#fig, ax = plt.subplots()
#fig.suptitle('kkk')
#k.showResults(ax=ax)
#print(os.path.dirname(os.path.realpath(__file__)))


#fresults = os.path.join('./test/f001_res.vtk')
#if os.path.isfile(fresults):
#    print('kk')
#    mesh_dict=mt.vtk_import(fresults)#makes a dictionary of a mesh 
#    mesh = Mesh_obj.mesh_dict2obj(mesh_dict)# this is a mesh_obj class instance 
#    mesh.show()
#


#%% test for timelapse inversion
#k = R2('/media/jkl/data/phd/tmp/r2gui/api/test/')
#k.createTimeLapseSurvey(os.path.join(k.dirname, 'testTimelapse'))
#k.invert(iplot=False)
#k.showSection(os.path.join(k.dirname, 'f001_res.vtk'))
#k.showSection(os.path.join(k.dirname, 'f002_res.vtk'))
