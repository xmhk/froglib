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



