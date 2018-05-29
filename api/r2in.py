import numpy as np
import os
import shutil

from mesh import Mesh

# all R2 parameters will be store in a dictionary to be written to R2.in file at the end of this script. All variables names are the same as in R2 manual.

class R2in:
    def __init__(self, param):
        ''' param is a dictionary that contains all the variables to be written to R2.in'''

    def writeR2(self, outputname=''):
        ''' write to R2.in
        '''

        if outputname == '':
            ouptutname = dirname + '/local' # write to local name