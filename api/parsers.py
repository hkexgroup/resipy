#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  1 11:23:23 2018

@author: jkl
"""

import numpy as np
import pandas as pd


def syscalParser(fname):
        df = pd.read_csv(fname, skipinitialspace=True)
        # delete space at the end and the beginning of columns names
        headers = df.columns
        newheaders = list(map(str.strip, headers)) 
        dico = dict(zip(headers, newheaders))
        df = df.rename(index=str, columns=dico)
        df = df.rename(columns={'Spa.1':'a',
                                'Spa.2':'b',
                                'Spa.3':'m',
                                'Spa.4':'n',
                                'In':'i',
                                'Vp':'vp',
                                'Dev.':'dev',
                                'M':'ip', #M1, M2,...Mn are good for now when importing
                                'Sp':'sp'})
        
        array = df[['a','b','m','n']].values
        spacing = np.unique(np.sort(array.flatten()))[1]
        array = np.round(array/spacing+1).astype(int)
        df[['a','b','m','n']] = array
        df['resist'] = df['vp']/df['i']
        imax = int(np.max(array))
        elec = np.zeros((imax,3))
        elec[:,0] = np.arange(0,imax)*spacing
                
        return elec, df
    
#%%test code
#elec, df = syscalParser('test/syscalFile.csv')

