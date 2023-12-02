import numpy as np
from resipy import Project

k = Project(typ='R3t')
elec = np.zeros((12, 3))
elec[:, 0] = np.arange(12)
elec[6:, 1] = 2
k.setElec(elec)
k.createMesh('tetra')
k.showMesh()
k.forward()
