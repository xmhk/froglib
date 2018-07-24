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


from matplotlib import pyplot as plt
import numpy as np
from .corefunctions import  frogtr
from .phasemanipulations import remove_shg_freq_ambi, remove_blind_freq_ambi

def compare_frog(mexp, pulse, gpulse, mode='shg',
                 fixshg_ambi=True,
                 fixblind_ambi=True):
    nn = len(pulse)
    mres = frogtr(pulse, gpulse, mode=mode)
    i1 = np.abs(mexp) ** 2 / np.sum(np.sum(np.abs(mexp) ** 2))
    i2 = np.abs(mres) ** 2 / np.sum(np.sum(np.abs(mres) ** 2))
    ferr = np.sqrt(1 / nn ** 2 * np.sum(np.sum(np.abs(i1 - i2) ** 2)))
    if mode=='shg' and fixshg_ambi:
        pulse = remove_shg_freq_ambi(pulse)
        gpulse = remove_shg_freq_ambi(gpulse)
    if mode=='blind' and fixblind_ambi:
        pulse, gpulse = remove_blind_freq_ambi(pulse, gpulse )

    plt.figure()
    plt.subplot(221)
    plt.imshow(i1)
    plt.xlabel("time")
    plt.ylabel("freq")
    plt.title("exp")
    plt.colorbar()
    plt.tight_layout()
    plt.subplot(222)
    ts = "rec, err = %.2e"%ferr
    plt.title(ts)
    plt.imshow(i2)
    plt.colorbar()
    plt.xlabel("time")
    plt.ylabel("freq")
    plt.tight_layout()
    plt.subplot(223)
    plt.imshow(i1 - i2)
    plt.xlabel("time")
    plt.ylabel("freq")
    plt.title("difference")
    plt.colorbar()
    plt.tight_layout()
    plt.subplot(224)
    plt.plot(np.abs(pulse), c='r', dashes=(3, 1), label='S')
    plt.plot(np.abs(gpulse), c='b', dashes=(3, 1), label='G')
    plt.legend(loc=0)
    plt.ylabel("amplitude")
    plt.xlabel("time")
    plt.twinx()
    plt.plot(np.unwrap(np.angle(pulse)), c='r')
    plt.plot(np.unwrap(np.angle(gpulse)), c='b')

    plt.ylabel("phase")
    plt.tight_layout()

def simplerecresult(M, rd):
    plt.figure(figsize=(5.5, 1.5))
    plt.semilogy(rd['errors'])
    plt.xlabel("iteration no")
    plt.ylabel("error")
    ts = "min error = %s"%rd['minerror']
    plt.title(ts)
    plt.tight_layout()
    compare_frog(M, rd['min_sp'], rd['min_gp'], mode=rd['mode'])