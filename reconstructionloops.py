from .corefunctions import pcgpstep
import numpy as np



def simplerec(mexp, iterations=10, mode='shg', pulse=None, gatepulse=None):
    errors = []
    sps = []
    gps = []
    for jj in range(iterations):
        pulse, gatepulse, frogerr = pcgpstep(mexp, pulse, gatepulse, mode=mode)
        errors.append(frogerr)
        sps.append(pulse)
        gps.append(gatepulse)
        # if not(extspec is None):    # use external spec as constraint? (not yet)
        #    nn = len(pulse); n2 = int(nn/2)
        #    PPS = np.roll(np.fft.fft(np.roll(pulse,n2)), n2)
        #    pulse = np.roll( np.fft.ifft(np.roll(np.abs(extspec[0]) * np.exp(1.0j * np.angle(PPS)),n2)), n2)
        #    GPS = np.roll(np.fft.fft(np.roll(gatepulse,n2)), n2)
        #    gatepulse = np.roll( np.fft.ifft( np.roll( np.abs(extspec[1]) * np.exp(1.0j * np.angle(GPS)), n2)), n2)
    rd = {'errors': errors, 'sp': sps, 'gp': gps, 'exp': mexp}
    # print(np.nonzero(np.isnan(rd['errors'])))
    minerr = np.min(rd['errors'])
    indx = np.nonzero(minerr == rd['errors'])[0][0]
    rd['minerror'] = minerr
    rd['min_sp'] = rd['sp'][indx]
    rd['min_gp'] = rd['gp'][indx]
    return rd
