#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 30 16:48:54 2018

@author: jkl
"""

import numpy as np
import pandas as pd
import os


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
    
    def readData(self, fname):
        ''' read electrode and quadrupoles data
        '''
        
    
    def parseData(self, fname):
        data = pd.read_csv(fname)
        return data
        
    def createR2in(self, param):
        ''' create R2.in object
        '''
        
    def invert(self):
        ''' perform inversion
        '''
    
    def show(self):
        ''' plot the results 
        '''
        
    def createMesh(self, typ='quad'):
        ''' create a mesh object
        typ:
            quad : quadrilateral mesh
            triang : triangular mesh
        '''
        
        self.mesh = mesh
        
        
    def showMesh(self):
        ''' call show mesh
        '''

        
#%% test code
k = R2()





       