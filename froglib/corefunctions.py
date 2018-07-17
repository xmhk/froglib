import numpy as np

def gaussianrandomphase(n):
    """Return a gaussian pulse with random phase.

    Arguments:

        n : length of desired gaussian pulse array

    Returns:

        field : gaussian pulse with random phase
                
    """
    tvec = np.arange(n) - int(n / 2)
    p = np.exp(-tvec ** 2 / (n / 16) ** 2) * np.exp(2.0j * np.pi * np.random.rand(n))
    #p = np.random.rand(n) * np.exp(2.0j * np.pi * np.random.rand(n))
    return p

def frogtr(field1, field2, mode='shg'):
    """Calculate the Frog Trace of given field(s).
    
    Arguments:

        field1, field2 : input fields

    Optional Arguments:

        mode : determines whether to calculate the SHG or the Blind Frog trace.

            'shg'=SHG trace
            'blind'=Blind Frog trace

    Returns:

         frogtrace : Field of the frog trace, which in general is complex.
                     To compare with some experimental (intensity) trace,
                     the absolute value squared of Frog trace has to be taken.
         
         
    """    
    nn = len(field1)
    n2 = int(nn / 2)
    if mode == 'shg':
        ap = np.outer(field1, field2) + np.outer(field2, field1)  # Axis: time-time
    elif mode == 'blind':
        ap = np.outer(field1, field2)
    m1 = np.zeros(np.shape(ap), dtype=np.complex128)
    m2 = np.zeros(np.shape(ap), dtype=np.complex128)
    for i in range(n2 - 1, -n2, -1):
        m1[i + n2, :] = np.roll(ap[i + n2, :], -i)
    m1 = np.transpose(m1)
    for i in range(nn):
        m2[i, :] = np.roll(np.fft.fft(np.roll(m1[i, :], +n2)), -n2)  # time-freq
    m2 = np.transpose(m2)  # freq - time
    m2 = m2 / np.max(np.max(np.abs(m2)))
    return m2


def pcgpstep(mexp, pulse, gatepulse, mode='shg', svd='full'):
    """Make a step in the pcgp algorithm.

    This function implements the PCGP Algoritm describe in Rick Trebinos Book

    Trebino, Rick. Frequency-resolved optical gating: the measurement of
    ultrashort laser pulses. Springer Science & Business Media, 2012.

    Arguments:

        mexp: amplitude represenation of the experimental FROG trace, e.g.
              the square root of the measured intensity.

        pulse: input pulse to use as signal. When pulse==None, a gaussian with
               random phase is used.

        gatepulse: input pulse to use as gate. When pulse==None, a gaussian with
                    random phase is used.

    Optional arguments:

        mode: can either be 'shg' or 'blind'. Determines the kind of Frog trace to calculate.

        svd: can either be 'full' or 'power'. When 'full' is used, the singular value decomposition
             of numpy is used. For 'power',  the (faster) 'power method' is used.

    Returns:

        pulse, gatepulse : iterated fields for signal and gate

        ferr : Frog error for the reconstructed trace

    """
    if pulse is None:
        nn = np.shape(mexp)[0]
        pulse = gaussianrandomphase(nn)
    if gatepulse is None:
        nn = np.shape(mexp)[0]
        gatepulse = gaussianrandomphase(nn)
    nn = len(pulse)
    n2 = int(nn / 2)
    m2 = frogtr(pulse, gatepulse, mode=mode)
    i1 = np.abs(mexp) ** 2 / np.sum(np.sum(np.abs(mexp) ** 2))
    i2 = np.abs(m2) ** 2 / np.sum(np.sum(np.abs(m2) ** 2))
    ferr = np.sqrt(1 / nn ** 2 * np.sum(np.sum(np.abs(i1 - i2) ** 2)))

    m3 = np.abs(mexp) * np.exp(1.0j * np.angle(m2))
    m3 = np.transpose(m3)  # zeit - freq
    m4 = np.zeros(np.shape(m2), dtype=np.complex128)
    m5 = np.zeros(np.shape(m2), dtype=np.complex128)

    for i in range(nn):
        m4[i, :] = np.roll(np.fft.ifft(np.roll(m3[i, :], -n2)), n2)
    for i in range(n2 - 1, -n2, -1):
        m5[i + n2, :] = np.roll(m4[:, i + n2], i)  # time-time
    if svd=='full':
        # full SVD
        u, w, v = np.linalg.svd(m5)
        pulse = u[:, 0]
        gatepulse = v[0, :]
    else:
        #  power method
        pulse = np.dot(np.dot(m5, np.transpose(m5)), pulse)
        gatepulse = np.dot(np.dot(np.transpose(m5), m5), gatepulse)
        pulse = pulse / np.sqrt( np.sum( np.abs(pulse)**2))
        gatepulse = gatepulse / np.sqrt(np.sum(np.abs(gatepulse) ** 2))
    return pulse, gatepulse, ferr