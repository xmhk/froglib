import os, sys
sys.path.append(os.path._getfullpathname("../"))

import numpy as np
from matplotlib import pyplot as plt
import matplotlib
from froglib import *

tvec = np.linspace(-8,8, 128)
#double pulse
t0=0.4
deltat = 0.7
signal1 = np.exp(-((tvec+deltat/2)/t0)**2) - np.exp(-((tvec-deltat/2)/t0)**2) * np.exp((-(tvec-deltat/2)/t0 * -1j))
signal2 = 1./np.cosh(tvec) * np.exp(1.0j * (0.6 * tvec + 0.5 * tvec**2))

# generate traces

F1 = frogtr(signal1, signal1, mode='shg')
F2 = frogtr(signal2, signal2, mode='shg')
F3 = frogtr(signal1, signal2, mode='blind')

# plot traces

spp = matplotlib.figure.SubplotParams(wspace=0.5)
plt.figure(figsize=(7,3), subplotpars=spp)
ax1= plt.subplot(131);plt.title("double gaussian, SHG")
plt.imshow(np.abs(F1)**2)
ax2= plt.subplot(132);plt.title("sech, SHG")
plt.imshow(np.abs(F2)**2)
ax3=plt.subplot(133);plt.title("double gaussian - sech,\nBLIND")
plt.imshow(np.abs(F3)**2)

for ax in [ax1,ax2,ax3]:
    ax.set_xlabel("time")
    ax.set_ylabel("freq")
plt.savefig("pics/example_traces.png")
plt.show()