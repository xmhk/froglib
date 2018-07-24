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


from .corefunctions import frogtr, pcgpstep
from .reconstructionloops import simplerec, mixfrog
from .statistics import deviationtau, calc_average
from .auxiliary import flip, shift_center

from .phasemanipulations import remove_relative_phase_offset
from .phasemanipulations import remove_shg_freq_ambi
from .phasemanipulations import remove_blind_freq_ambi
from .phasemanipulations import specshift_to_ref

from .plotting import compare_frog, simplerecresult
