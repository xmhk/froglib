import numpy as np

def remove_relative_phase_offset(field1, field2, border=6,
                          intermediate_steps=50,
                           magnitudenitudes=np.linspace(0, 4, 5)):
    """remove the arbitrary phase differences between two fields (linar+offset)"""
    nn = len(field1)
    n2 = int(nn / 2.0)
    xv = np.arange(nn) - n2
    testSP1 = np.roll(np.fft.fft(np.roll(field1, n2)), n2)
    testSP2 = np.roll(np.fft.fft(np.roll(field2, n2)), n2)
    NN = intermediate_steps
    MM = np.zeros((NN, NN))
    PosNullX = 0
    PosNullY = 0
    for magnitude in magnitudenitudes:
        jjv = np.linspace(-border * 10 ** -magnitude, border * 10 ** -magnitude, NN)\
              + PosNullX
        kkv = np.linspace(-border * 10 ** -magnitude, border * 10 ** -magnitude, NN)\
              + PosNullY
        for j, jj in enumerate(jjv):
            for k, kk in enumerate(kkv):
                MM[j, k] = np.sum(
                    (np.imag(testSP1) - np.imag(testSP2 \
                                * np.exp(1.0j * (xv * jj + kk)))) ** 2
                    * (np.real(testSP1) - np.real(testSP2 \
                                * np.exp(1.0j * (xv * jj + kk)))) ** 2
                )
        minindex = np.nonzero(MM == np.min(np.min(MM)))
        PosNullX = jjv[minindex[0][0]]
        PosNullY = kkv[minindex[1][0]]

    field3 = np.roll(
        np.fft.ifft(
            np.roll(
                testSP2 * np.exp(1.0j * (PosNullX * xv + PosNullY))
                , n2)
        ), n2)
    return field3


def remove_shg_freq_ambi(pulse):
    """test and remove the shg frequency ambiguity (spek shift by n/2)"""
    l = len(pulse)
    l2 = int(l / 2)
    tests = np.exp(- ((np.arange(l) - l2) / (l / 4)) ** 6)
    pulsspeksh = np.fft.ifft(np.fft.fftshift(np.fft.fft(pulse)))
    normal = np.sum(np.abs(np.multiply(tests, np.fft.fftshift(np.fft.fft(pulse)))))
    shifted = np.sum(np.abs(np.multiply(tests, np.fft.fftshift(np.fft.fft(pulsspeksh)))))
    if normal > shifted:
        return pulse
    else:
        return pulsspeksh



def remove_blind_freq_ambi(pulse, gpulse, anz=200):
    """experimental"""
    errors = []
    for i in range(anz):
        ll = len(pulse)
        xv = np.arange(ll) - int(ll / 2)
        dom = np.pi / anz * 2 * i
        sps = np.fft.fftshift(np.abs(np.fft.fft(pulse * np.exp(1.0j * dom * xv))))
        spg = np.fft.fftshift(np.abs(np.fft.fft(gpulse * np.exp(-1.0j * dom * xv))))
        errors.append(np.sum(np.abs(sps - spg) ** 2))
    minerr = np.min(errors)
    ii = np.nonzero(minerr == errors)[0][0]
    dom = np.pi / anz * 2 * ii
    pulse = pulse * np.exp(1.0j * dom * xv)
    gpulse = gpulse * np.exp(-1.0j * dom * xv)
    return pulse, gpulse



def specshift_to_ref(referencepulse, pulse, anz=200):
    """shift spectrum of some pulse for maximum overlap with referencepuse spec"""
    errors = []
    for i in range(anz):
        ll = len(referencepulse)
        xv = np.arange(ll) - int(ll / 2)
        dom = np.pi / anz * 2 * i
        sps = np.fft.fftshift(np.abs(np.fft.fft(referencepulse)))
        spg = np.fft.fftshift(np.abs(np.fft.fft(pulse * np.exp(-1.0j * dom * xv))))
        errors.append(np.sum(np.abs(sps - spg) ** 2))
    minerr = np.min(errors)
    ii = np.nonzero(minerr == errors)[0][0]
    dom = np.pi / anz * 2 * ii
    pulse = pulse * np.exp(-1.0j * dom * xv)
    return pulse