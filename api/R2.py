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
from subprocess import PIPE, call, Popen
OS = platform.system()

sys.path.append(os.path.relpath('..'))

#from QuadMesh import QuadMesh
#from TriMesh import TriMesh

from api.Survey import Survey
from api.r2in import write2in

class R2(object): # R2 master class instanciated by the GUI
    def __init__(self, dirname=''):
        if dirname == '':
            dirname = os.getcwd()
            print('using the current directory:', dirname)
        self.dirname = dirname # working directory
        self.surveys = [] # list of survey object
        self.surveysInfo = [] # info about surveys (date)
        self.mesh = None # mesh object (one per R2 instance)
        
        
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
        
    
    def createMesh(self, typ='quad'):
        ''' create a mesh object
        typ:
            quad : quadrilateral mesh
            triang : triangular mesh
        '''
        if typ == 'quad':
            mesh = QuadMesh(elec, nnode=4)
        if typ == 'trian':
            mesh = TriMesh(elec, resolution=1)
        self.mesh = mesh
        
        
    def invert(self, param={}):
        ''' invert the data, first generate R2.in file, then run
        inversion using appropriate wrapper, then return results
        '''
        # write configuration file
#        write2in(param, self.dirname)
        
        # copy R2.exe
#        os.copy('../external-exe/R2.exe',self.dirname)

        self.runR2()
        
    def runR2(self):
        # run R2.exe
        cwd = os.getcwd()
        os.chdir(self.dirname)
        
        if OS == 'Linux':
            cmd = ['wine','R2.exe']
        if OS == 'Windows':
            cmd = ['R2.exe']
        
        p = Popen(cmd, stdout=PIPE, shell=False)
        while p.poll() is None:
            line = p.stdout.readline().rstrip()
            print(line.decode('utf-8'))
            
        os.chdir(cwd)
        
        
#%% test code
#k = R2('/media/jkl/data/phd/tmp/r2gui/api/test')
#k.createSurvey('test/syscalFile.csv', ftype='Syscal')
#k.pseudo(contour=True)
#k.linfit(iplot=True)
#k.invert()
