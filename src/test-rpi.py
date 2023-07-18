import numpy as np
from resipy import Project

k = Project()
elec = np.zeros((12, 3))
elec[:, 0] = np.arange(12)
k.setElec(elec)
k.createMesh('trian')
k.forward()
