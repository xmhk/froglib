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

# two gaussians of different height and width, different freqs
# chirped gaussian as gate pulse

field = 0.5 * np.exp( - (tvec-deltat)**2) - np.exp( - (tvec+deltat)**2/0.5**2) * np.exp(1.0j * freqoff * tvec)
gate= np.exp(-tvec**2/0.16) * np.exp(-1.3j * tvec**2)

Mexp = np.abs(frogtr(field, gate, mode='blind'))


res = simplerec(Mexp, iterations=50, mode='blind')

# show result

simplerecresult(Mexp, res)

plt.show()



