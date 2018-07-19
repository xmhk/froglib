from .corefunctions import pcgpstep
import numpy as np
from matplotlib import pyplot as plt
from .plotting import simplerecresult



def simplerec(mexp, iterations=10, mode='shg', pulse=None, gatepulse=None, svd='full'):
    """Simple reconstruction loop for Frog traces.

    Arguments:

        mexp : amplitude represenation of the experimental FROG trace, e.g.
               the square root of the measured intensity.

    Optional Arguments:

        iterations : number of PCGP steps to iterate (default = 10)

        mode : may be 'shg' or 'blind'  (default = 'shg')

        pulse : when given, this pulse will be feed into the PCGP loop (default = None)

        gatepulse : when given, this gatepulse will be feed into the PCGP loop (default = None)

        svd : method to calculate the singular value decomposition. 'full' : use numpy's SVD,
              'power' : use the power method.

    Returns:

        rdict : dictionary holding the values:

            errors : Array holding Frog error for each iteration

            gp, sp :  Arrays with signal and gate fields for each iteration

            minerror : minimal Frog error that occured during iterations

            min_sp, min_gp : signal and gate pulses for minimal Frog error

            mode : mode used ('shg' or 'blind')

            exp : experimental trace (amplitude)

    """

    errors = []
    sps = []
    gps = []
    for jj in range(iterations):
        pulse, gatepulse, frogerr = pcgpstep(mexp, pulse, gatepulse, mode=mode, svd=svd)
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
    minerr = np.min(rd['errors'])
    indx = np.nonzero(minerr == rd['errors'])[0][0]
    rd['minerror'] = minerr
    rd['min_sp'] = rd['sp'][indx]
    rd['min_gp'] = rd['gp'][indx]
    rd['mode']=mode
    return rd



def mixfrog(mexp, startnum=3, startite=20, itestep2=100, plot=False):
    """Try different reconstructions, return the best.
    
    This algorithm makes a number of different attempts to reconstruct 
    a trace, both with 'shg' and 'blind' mode.
    
    The best result (minimal Frog error) is iterated further.
    
    Arguments: 

        mexp : experimental Frog trace (amplitude)

    Optional arguments:

        startnum : number of attempts for each 'shg' and 'blind' mode

        startite : number of iterations to perform for each mode

        itestep2 : number of iterations for the best field

        plot : whether to plot some results for this mode.

    Returns:

        rdict : dictionary holding the values:

            errors : Array holding Frog error for each iteration

            gp, sp :  Arrays with signal and gate fields for each iteration

            minerror : minimal Frog error that occured during iterations

            min_sp, min_gp : signal and gate pulses for minimal Frog error

            mode : mode used ('shg' or 'blind')

            exp : experimental trace (amplitude)

    """
    
    errors = []
    pp = []
    gp = []
    for i in range(startnum):
        rd = simplerec(mexp, iterations=startite, mode='shg')
        errors.append(rd['minerror'])
        pp.append(rd['min_sp'])
        gp.append(rd['min_gp'])
    for i in range(startnum):
        rd = simplerec(mexp, iterations=startite, mode='blind')
        errors.append(rd['minerror'])
        pp.append(rd['min_sp'])
        gp.append(rd['min_gp'])

    minerr = np.min(errors)
    iindx = np.nonzero(np.array(errors) == np.min(errors))[0][0]
    if iindx < startnum:
        mode = 'shg'
    else:
        mode = 'blind'

    rd = simplerec(mexp, pulse=pp[iindx], gatepulse=gp[iindx], mode=mode,
                   iterations=itestep2)
    if plot:
        plt.figure(figsize=(5,2))
        plt.plot(errors, 'o')
        plt.axvline(x=startnum-1, lw=0.8, c='0.0', dashes=(5,1))
        plt.xlabel("Round no")
        plt.ylabel("Frog error")
        simplerecresult(mexp, rd)
    return rd

