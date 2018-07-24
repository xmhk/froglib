# froglib

A library to reconstruct the complex optical field of short optical pulses from FROG 
(Frequency Resolved Optical Gating[1]) traces.

For a comprehensive overview on the FROG technique, see: 
```
  Trebino, Rick. Frequency-resolved optical gating: the measurement of ultrashort laser pulses. 
  Springer Science & Business Media, 2012.
```
### What does froglib provide?

 - generate artificial FROG traces for Second Harmonic Generation (SHG) and Blind FROG setup
 - basic implementation of the PCGP (Principle Components Generalized
Projection) algorithm
 - basic reconstruction loops 
 - basic functions to remove simple ambiguities
 - tools to calculate the field with errorbars from a list of reconstructed fields
 
 ### What DOESN'T it provide?
  
 - A full-feature, ready-to-use end user program. In general, several unforeseeable things may happen
   when reconstructing a FROG trace. Even simple ambiguities are sometimes hard to detect automatically.
   Reconstruction sometimes tends to stick at traces with small Frog errors, which for the human inspector
    are clearly different from the experimental traces. In general, human interaction and judgement is unavoidable. 

 - Automatic calibration. Time and frequency grid of a FROG trace have to be choosen according to the Fourier 
   relation. Possibly, measure traces have to be stretched accordingly before starting reconstruction.  
   All calculations in froglib are done in normalized units (pixel).
   

## Usage

### generate FROG traces: **frogtr**
 ![example frog traces](examples/pics/example_traces.png)    
```
    import numpy as np
    from matplotlib import pyplot as plt
    import matplotlib
    from froglib import *
    
    tvec = np.linspace(-8,8, 128)
    #double pulse
    t0=0.4
    deltat = 0.7
    signal1 = np.exp(-((tvec+deltat/2)/t0)**2) - np.exp(-((tvec-deltat/2)/t0)**2) * np.exp((-(tvec-deltat/2)/t0 * -1j))
    signal2 = 1./np.cosh(tvec) * np.exp(1.0j * (0.6 * tvec + 0.5 * tvec**2))
    
    # generate traces
    
    F1 = frogtr(signal1, signal1, mode='shg')
    F2 = frogtr(signal2, signal2, mode='shg')
    F3 = frogtr(signal1, signal2, mode='blind')
    
    # plot traces
    
    spp = matplotlib.figure.SubplotParams(wspace=0.5)
    plt.figure(figsize=(7,3), subplotpars=spp)
    ax1= plt.subplot(131);plt.title("double gaussian, SHG")
    plt.imshow(np.abs(F1)**2)
    ax2= plt.subplot(132);plt.title("sech, SHG")
    plt.imshow(np.abs(F2)**2)
    ax3=plt.subplot(133);plt.title("double gaussian - sech,\nBLIND")
    plt.imshow(np.abs(F3)**2)
    
    for ax in [ax1,ax2,ax3]:
        ax.set_xlabel("time")
        ax.set_ylabel("freq")
```
   



   

### simple, single reconstruction: **simplerec**
 ![example reconstruction error](examples/pics/rec_artif_shg0.png)
 ![example reconstruction](examples/pics/rec_artif_shg.png)
```
import numpy as np
from matplotlib import pyplot as plt
from froglib import *

# generate signal
tmax = 7
N = 128
tvec = np.linspace(-tmax, tmax, N)
deltat = 0.5
freqoff = 1.2
#two gaussians of different height and width, different freqs
field = 0.5 * np.exp( - (tvec-deltat)**2) - np.exp( - (tvec+deltat)**2/0.5**2) * np.exp(1.0j * freqoff * tvec)

# generate artifical trace (SHG)
Mexp = np.abs(frogtr(field, field, mode='shg'))

# use a simple loop of the PCGP algorithm
res = simplerec(Mexp, iterations=150)

# show result
simplerecresult(Mexp, res)
plt.show()
```

The result of the reconstruction is stored in the dictionary *res*:

``` 
    (keys) : (values)
    errors : Array holding Frog error for each iteration
    gp, sp :  Arrays with signal and gate fields for each iteration
    minerror : minimal Frog error that occured during iterations
    min_sp, min_gp : signal and gate pulses for minimal Frog error
    mode : mode used ('shg' or 'blind')
    exp : experimental trace (amplitude)
```


### more advanced reconstruction loop: mixfrog

The basic idea behind this is to make several attempts to reconstruct a trace.
Each attempt initially is only iterated some steps. Then, the best 'candidate' is 
chosen and iterated fully.
Both SHG and Blind FROG reconstructions are attempted. Usually, due to the fact that
experimental SHG FROG traces show a slight asymmetric shape, be it due to noise or misalignment,
the Blind part will give better result.

When SHG traces are reconstructed with in 'blind' mode, both reconstructed signal and gate pulse
will look very similar albeit phase/time ambiguities.

For more information, see 

```
Hause, A., et al. "Reliable multiple-pulse reconstruction from second-harmonic-generation 
frequency-resolved optical gating spectrograms." JOSA B 32.5 (2015): 868-877.
``` 

![error for different attempts](examples/pics/mix0.png)
![error for different attempts](examples/pics/mix1.png)
![error for different attempts](examples/pics/mix2.png)
```
import numpy as np
from matplotlib import pyplot as plt
from froglib import *


M = np.loadtxt("example_data/stronglychirpeddouble.dat")
Mexp = np.sqrt(M)
res = mixfrog(Mexp, startnum=10, plot=True )

plt.show()
```


## License

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

* Christoph Mahnke, 2018      (chmhnk_at_googlemail_dot_com)