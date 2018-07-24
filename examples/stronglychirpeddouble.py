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


M = np.loadtxt("example_data/stronglychirpeddouble.dat")
Mexp = np.sqrt(M)
res = mixfrog(Mexp, startnum=10, plot=True )

plt.show()



