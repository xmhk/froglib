import numpy as np
from matplotlib import pyplot as plt
from .corefunctions import frogtr, pcgpstep
from .reconstructionloops import rek1



def shiftmitte(feld, eps05=5e-2):
    pp = np.abs(feld) ** 2
    cs = np.cumsum(pp)
    rel = cs / np.sum(pp)
    aa = np.nonzero(np.multiply(rel > 0.5 - eps05, rel < 0.5 + eps05))[0]
    ll = len(feld)
    shift = 0
    if len(aa) > 0:
        shift = int(ll / 2 - np.round(np.mean(aa)))
    return np.roll(feld, shift)


def shg_check_freq_ambiguity(pulse):
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


def blindshift(pulse, gpulse, anz=200):
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
    # print(ii)
    # print(errors)
    # plt.figure()
    # plt.plot(np.fft.fftshift(np.abs(np.fft.fft(pulse))))
    # plt.plot(np.fft.fftshift(np.abs(np.fft.fft(gpulse))))
    pulse = pulse * np.exp(1.0j * dom * xv)
    gpulse = gpulse * np.exp(-1.0j * dom * xv)
    pulse = shg_check_freq_ambiguity(pulse)
    gpulse = shg_check_freq_ambiguity(gpulse)
    # plt.figure()
    # plt.plot(np.fft.fftshift(np.abs(np.fft.fft(pulse))))
    # plt.plot(np.fft.fftshift(np.abs(np.fft.fft(gpulse))))
    return pulse, gpulse





def frogvergleich(mexp, pulse, gpulse, mode='shg'):
    nn = len(pulse)
    mres = frogtr(pulse, gpulse, mode=mode)
    i1 = np.abs(mexp) ** 2 / np.sum(np.sum(np.abs(mexp) ** 2))
    i2 = np.abs(mres) ** 2 / np.sum(np.sum(np.abs(mres) ** 2))
    ferr = np.sqrt(1 / nn ** 2 * np.sum(np.sum(np.abs(i1 - i2) ** 2)))
    plt.figure()
    plt.subplot(221)
    plt.imshow(i1)
    plt.colorbar()
    plt.subplot(222)
    plt.imshow(i2)
    plt.colorbar()
    plt.subplot(223)
    plt.imshow(i1 - i2)
    plt.colorbar()
    plt.subplot(224)
    plt.plot(np.abs(pulse), c='r', dashes=(3, 1))
    plt.plot(np.abs(gpulse), c='b', dashes=(3, 1))
    plt.twinx()
    plt.plot(np.unwrap(np.angle(pulse)), c='r')
    plt.plot(np.unwrap(np.angle(gpulse)), c='b')
    print(ferr)
