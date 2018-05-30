#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 30 17:02:37 2018

@author: jkl
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt # the only script allowed to import matplotlib

from r2 import R2


class R2m(R2):
    def __init__(self, dirname=''):
        R2.__init__(self, dirname)
    
    def show(self):
        ''' plot section
        '''
        fig, ax = plt.subplots()
        ax.imshow(np.random.randn(100,40))
        fig.show()


#%% test code

k = R2m()
k.show()

