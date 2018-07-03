import numpy as np

def flip(v):
    """flips a vector v"""
    v2 = [x for x in v]
    v2.reverse()
    return np.array(v2)



def shift_center(feld, eps05=5e-2):
    """shift some field to the center of the time window, with respect to energy"""
    pp = np.abs(feld) ** 2
    cs = np.cumsum(pp)
    rel = cs / np.sum(pp)
    aa = np.nonzero(np.multiply(rel > 0.5 - eps05, rel < 0.5 + eps05))[0]
    ll = len(feld)
    shift = 0
    if len(aa) > 0:
        shift = int(ll / 2 - np.round(np.mean(aa)))
    return np.roll(feld, shift)