#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 30 16:48:54 2018

@author: jkl
"""

import numpy as np
import pandas as pd
import os

from QuadMesh import QuadMesh
from TriMesh import TriMesh


class R2():
    def __init__(self, dirname=''):
        if dirname == '':
            dirname = os.getcwd()
            print('using the current directory:', dirname)
        self.dirname = dirname
        
        
    def setwd(self, dirname):
        ''' set the working directory
        '''
        self.dirname = self.dirname
    
    
    def creatSurvey(self, fname):
        ''' read electrode and quadrupoles data and return 
        a survey object
        '''    
    
    
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
        
        
    def invert(self, param):
        ''' invert the data, first generate R2.in file, then run
        inversion using appropriate wrapper, then return results
        '''
        
        
        
        
#%% test code
k = R2()





       