import numpy as np
from matplotlib import pyplot as plt

from .corefunctions import frogtr, pcgpstep
from .reconstructionloops import simplerec
from .statistics import deviationtau
from .auxiliary import flip, shift_center

from .phasemanipulations import remove_relative_phase_offset
from .phasemanipulations import remove_shg_freq_ambi
from .phasemanipulations import remove_blind_freq_ambi
from .phasemanipulations import specshift_to_ref

from .plotting import compare_frog, simplerecresult

#
# scratch functions -> sort or replace
#


""" 
def mixtest(M, startnum=5, startite=30):
    errors = []
    pp = []
    gp = []
    for i in range(startnum):
        rd = simplerec(M, iterations=startite, mode='shg')
        errors.append(rd['minerror'])
        pp.append(rd['min_sp'])
        gp.append(rd['min_gp'])
    for i in range(startnum):
        rd = simplerec(M, iterations=startite, mode='blind')
        errors.append(rd['minerror'])
        pp.append(rd['min_sp'])
        gp.append(rd['min_gp'])
    plt.figure()
    plt.plot(errors)
    minerr = np.min(errors)
    iindx = np.nonzero(np.array(errors)==np.min(errors))[0][0]
    rd = simplerec(M, pulse=pp[iindx], gatepulse=gp[iindx], iterations=100)
    rd2 = simplerec(M, pulse=pp[iindx], gatepulse=gp[iindx], mode='blind', iterations=100)
    frogvergleich(M, rd['min_sp'], rd['min_gp'])
    frogvergleich(M, rd2['min_sp'], rd2['min_gp'], mode='blind')
"""











"""
