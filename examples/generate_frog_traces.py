"""
This file is part of froglib, which allows the creation
and reconstruction of FROG traces.

froglib is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

froglib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with froglib.  If not, see <http://www.gnu.org/licenses/>.

Contributors:

Christoph Mahnke, 2018

"""

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
plt.show()