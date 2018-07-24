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


import numpy as np

def flip(v):
    """Flips (reverses) a vector v.

    Arguments:

        v : vector

    Returns:

        v2 : reversed vector
        
    """
    v2 = [x for x in v]
    v2.reverse()
    return np.array(v2)



def shift_center(field, eps05=5e-2):
    """Shift some field to the center of the time window, with respect to energy.
    
    This function shifts (rotates) some field to be centered in the time window,
    with respect to energy.
    
    Arguments:
          
        field : input field

    Optional arguments:

        eps05 : defines how the center of the field is determined.
                The function search for the field array element which
                fulfills (cumulated energy) = | 0.5 - eps05 | x total energy.

    Returns:

        shiftfield : shifted version of the field

    """
    pp = np.abs(field) ** 2
    cs = np.cumsum(pp)
    rel = cs / np.sum(pp)
    aa = np.nonzero(np.multiply(rel > 0.5 - eps05, rel < 0.5 + eps05))[0]
    ll = len(field)
    shift = 0
    if len(aa) > 0:
        shift = int(ll / 2 - np.round(np.mean(aa)))
    return np.roll(field, shift)