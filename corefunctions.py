import numpy as np

def gaussianrandomphase(n):
    tvec = np.arange(n) - int(n / 2)
    p = np.exp(-tvec ** 2 / (n / 16) ** 2) * np.exp(2.0j * np.pi * np.random.rand(n))
    #p = np.random.rand(n) * np.exp(2.0j * np.pi * np.random.rand(n))
    return p

def frogtr(feld1, feld2, mode='shg'):
    nn = len(feld1)
    n2 = int(nn / 2)
    if mode == 'shg':
        ap = np.outer(feld1, feld2) + np.outer(feld2, feld1)  # Axis: time-time
    elif mode == 'blind':
        ap = np.outer(feld1, feld2)
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


def pcgpstep(mexp, pulse, gatepulse, mode='shg'):
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
    if True:
        # full SVD
        u, w, v = np.linalg.svd(m5)
        pulse = u[:, 0]
        gatepulse = v[0, :]
    else:
        #  power method
        pulse = np.dot(np.dot(m5, np.transpose(m5)), pulse)
        gatepulse = np.dot(np.dot(np.transpose(m5), m5), gatepulse)
    return pulse, gatepulse, ferr