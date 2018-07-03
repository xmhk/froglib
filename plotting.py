
from matplotlib import pyplot as plt
import numpy as np
from .corefunctions import  frogtr

def compare_frog(mexp, pulse, gpulse, mode='shg'):
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
    ts = "Ferr = %.2e"%ferr
    plt.title(ts)
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


def simplerecresult(M, rd):
    plt.figure(figsize=(5.5, 1.5))
    plt.semilogy(rd['errors'])
    ts = "min = %s"%rd['minerror']
    plt.title(ts)
    plt.tight_layout()
    compare_frog(M, rd['min_sp'], rd['min_gp'], mode=rd['mode'])