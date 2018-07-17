import numpy as np

def deviationtau(n):
    """Return the prefactor for error margins of n-point sample.

    Depending on the length of some series of measurements, the error
    margin is calculated by taking the standard deviation times a prefactor
    tau(n) / sqrt (n).

    error = tau(n) / sqrt (n) x standard_deviation

    This function calculates tau(n) / sqrt (n) and returns it.

    Arguments:

         n : number of samples

    Returns:

        prefactor

    """
    if n<=2:
        n = 2
    P = [12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262,
         2.228, 2.201, 2.179, 2.160, 2.145, 2.131, 2.120, 2.110, 2.101,
         2.093, 2.086, 2.080, 2.074, 2.069, 2.064, 2.060, 2.056, 2.052,
         2.048, 2.045, 2.042]
    if n > 31:
        return 2.0 / np.sqrt(n)
    else:
        return P[n - 2] / np.sqrt(n)


def calc_average(fieldlist):
    """Calculates the averages and errors for a list of fields.

    Arguments:

        fieldlist: list of input fields

    Returns:

        rdict: dictionary containing the values:

            field_mean : mean field (complex)

            re_mean : mean of the real part

            re_error : error for the real part

            im_mean : mean of the imag part

            im_error : error for the imag part

            amp_mean : mean amplitude of the field (real)

            amp_error : error for the amp_mean

            phase_mean : mean of the phase

            phase_error : error for the phase.

    """
    Mfieldlist = np.array(fieldlist)
    n2 = int(len(fieldlist[0]) / 2)
    # real part
    mre = np.mean(np.real(Mfieldlist), axis=0)
    stdre = np.std(np.real(Mfieldlist), axis=0)
    errorre = deviationtau(len(fieldlist)) * stdre
    # imag part
    mim = np.mean(np.imag(Mfieldlist), axis=0)
    stdim = np.std(np.imag(Mfieldlist), axis=0)
    errorim = deviationtau(len(fieldlist)) * stdim
    # abs
    meanfield = mre + 1.0j * mim
    # phase
    phase = np.unwrap(np.angle(meanfield))
    phase = phase - phase[n2]
    errorabs = np.abs(mre / np.sqrt(mre ** 2 + mim ** 2) * errorre) \
                + np.abs(mim / np.sqrt(mre ** 2 + mim ** 2) * errorim)
    errorphase = np.abs(mim / ((1.0 + mim ** 2 / mre ** 2) * mre ** 2) * errorre) \
                  + np.abs(mim / ((1.0 + mim ** 2 / mre ** 2) * mre ** 2) * errorim)
    rd = {}
    rd['field_mean'] = meanfield
    rd['re_mean'] = mre
    rd['re_error'] = errorre
    rd['im_mean'] = mim
    rd['im_error'] = errorim
    rd['amp_mean'] = np.abs(meanfield)
    rd['amp_error'] = errorabs
    rd['phase_mean'] = phase
    rd['phase_error'] = errorphase

    return rd