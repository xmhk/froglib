"""
This file is part of froglib, which allows the creation
and reconstruction of FROG traces.

froglib is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

froglib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with froglib.  If not, see <http://www.gnu.org/licenses/>.

Contributors:

Christoph Mahnke, 2018

"""

import numpy as np

def remove_relative_phase_offset(field1, field2, border=6,
                          intermediate_steps=50,
                           magnitudes=np.linspace(0, 4, 5)):
    """Remove the arbitrary phase differences between two fields (linear+offset).

    The relative phases for two fields in the context of Frog traces can differ by some
    arbitrary offsets: the first one is the constant offset, the second one is the slope
    of the phase in the spectral domain (which corresponds to the temporal position).
    This algorithm removes these two for the second input field with respect to the first one.
    It varies the spectral phase (offset + linear component) and tries to find the
    best overlap.

    Arguments:

        field1 : reference field

        field2 : field whose phase shall be adjusted;

    Optional arguments:

        border : Numerical value for the phases to be varied in between

        intermediate_steps : number of intermediate steps for each magnitude

        magnitudes: np.array or list of powers of 10. Phase offsets and slopes
                    as border * 10 **(-magnitudes) are tried and refines successively.

    Returns:

        field3 : field2 with adjusted phase.

    """
    nn = len(field1)
    n2 = int(nn / 2.0)
    xv = np.arange(nn) - n2
    testSP1 = np.roll(np.fft.fft(np.roll(field1, n2)), n2)
    testSP2 = np.roll(np.fft.fft(np.roll(field2, n2)), n2)
    NN = intermediate_steps
    MM = np.zeros((NN, NN))
    PosNullX = 0
    PosNullY = 0
    for magnitude in magnitudes:
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
    """Test and remove the shg frequency ambiguity (spectral shift by n/2).

    During reconstruction of SHG Frog, a common ambiguity occuring is a field shifted by
    half of the window in the spectral domain. As it gives the same (often correct) Frog
    trace as the unshifted field, a physically valid field can be constructed by shifting
    the spectrum back to the center of the frequency domain.

    Arguments:

         pulse : pulse field

    Returns:

        pulse2 : shifted pulse field

    """
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



def remove_blind_freq_ambi(pulse, gpulse, inter_steps=200, check_shg_ambi=True):
    """Experimental: remove blind frog frequency ambiguity."""
    errors = []
    for i in range(inter_steps):
        ll = len(pulse)
        xv = np.arange(ll) - int(ll / 2)
        dom = np.pi / inter_steps * 2 * i
        sps = np.fft.fftshift(np.abs(np.fft.fft(pulse * np.exp(1.0j * dom * xv))))
        spg = np.fft.fftshift(np.abs(np.fft.fft(gpulse * np.exp(-1.0j * dom * xv))))
        errors.append(np.sum(np.abs(sps - spg) ** 2))
    minerr = np.min(errors)
    ii = np.nonzero(minerr == errors)[0][0]
    dom = np.pi / inter_steps * 2 * ii
    pulse = pulse * np.exp(1.0j * dom * xv)
    gpulse = gpulse * np.exp(-1.0j * dom * xv)
    if check_shg_ambi:
        pulse = remove_shg_freq_ambi(pulse)
        gpulse = remove_shg_freq_ambi(gpulse)
    return pulse, gpulse



def specshift_to_ref(referencepulse, pulse, inter_steps=200):
    """Shift spectrum of some pulse for maximum overlap with reference pulse spectrum.
    
    Input:
    
        referencepulse : reference
        
        pulse : pulse whose spectrum shall be shifted
        
    Optional Arguments: 
    
        inter_steps : number of intermediate steps to try 
    
    """
    errors = []
    for i in range(inter_steps):
        ll = len(referencepulse)
        xv = np.arange(ll) - int(ll / 2)
        dom = np.pi / inter_steps * 2 * i
        sps = np.fft.fftshift(np.abs(np.fft.fft(referencepulse)))
        spg = np.fft.fftshift(np.abs(np.fft.fft(pulse * np.exp(-1.0j * dom * xv))))
        errors.append(np.sum(np.abs(sps - spg) ** 2))
    minerr = np.min(errors)
    ii = np.nonzero(minerr == errors)[0][0]
    dom = np.pi / inter_steps * 2 * ii
    pulse = pulse * np.exp(-1.0j * dom * xv)
    return pulse