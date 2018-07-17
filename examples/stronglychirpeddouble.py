import os, sys
sys.path.append(os.path._getfullpathname("../"))

import numpy as np
from matplotlib import pyplot as plt
from froglib import *


M = np.loadtxt("stronglychirpeddouble.dat")
Mexp = np.sqrt(M)

#res = simplerec(Mexp, iterations=50, mode='blind')
res = mixfrog(Mexp, startnum=10, plot=True )

# show result

#simplerecresult(Mexp, res)

plt.show()



