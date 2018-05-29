import numpy as np
import os

''' create a Mesh class for handling meshes

'''


class Mesh:
    def __init__(self, nbElectrode, spacing):
        ''' nbElectrode = number of electrodes
            spacing = spacing between two electrodes [m]
        '''
        self.nbElectrode = nbElectrode
        self.spacing = spacing
    
    def createRectangularMesh(self, nodesBetween=4):
        '''
        nodesBetween = number of notes between two electrodes
        '''


    def createTriangularMesh(self, geoInfo):
    ''' geoInfo = dict
    '''


    def write2vtk(self):
    ''' write a vtk version of the mesh
    '''

    

