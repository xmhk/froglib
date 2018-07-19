import os, sys
sys.path.append(os.path._getfullpathname("../"))

import numpy as np
from matplotlib import pyplot as plt
from froglib import *


M = np.loadtxt("example_data/stronglychirpeddouble.dat")
Mexp = np.sqrt(M)
res = mixfrog(Mexp, startnum=10, plot=True )

plt.show()



