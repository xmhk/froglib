import os, sys
sys.path.append(os.path._getfullpathname("../"))

import numpy as np
from matplotlib import pyplot as plt
from froglib import *


tmax = 7
N = 128

tvec = np.linspace(-tmax, tmax, N)

deltat = 0.5
freqoff = 1.2

#two gaussians of different height and width, different freqs

field = 0.5 * np.exp( - (tvec-deltat)**2) - np.exp( - (tvec+deltat)**2/0.5**2) * np.exp(1.0j * freqoff * tvec)


Mexp = np.abs(frogtr(field, field, mode='shg'))


res = simplerec(Mexp, iterations=50)

# show result

simplerecresult(Mexp, res)

plt.show()



