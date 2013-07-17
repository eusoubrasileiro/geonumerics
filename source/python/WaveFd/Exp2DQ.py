r"""

4th order explicit space test question

Convergence
 
.. math::
 
     \Delta t \leq \frac{2 \Delta s}{ V \sqrt{\sum_{a=-N}^{N} (|w_a^1| + |w_a^2|)}}

Where w_a are the centered differences weights

Under limits.

"""

import numpy as np
import sys

class RickerSource(object):
    r"""
    Ricker Wavelet:
    :math:`A = (1-2 \pi^2 f^2 t^2) e^{-\pi^2 f^2 t^2}`

    * fc        : maximum desired frequency
    * amp       : source strength 
    * delay     : to turn the source causal recommend 1/fc (starting at zero)
    """

    def __init__(self, fc, amp, delay=0.0):
        self.fc = fc
        self.amp = amp
        self.delay = delay

    def __call__(self, t):
        t = t-self.delay
        cx = (np.pi*self.fc*t)**2
        return  self.amp*(1-2*cx)*np.exp(-cx) 

Nz = Nx = 80
Dt = 0.002
Ds = 8.0
numberiter = 300


Uprevious = np.zeros([Nz, Nx]) # previous time
Ucurrent = np.zeros([Nz, Nx]) # current time
Ufuture = np.zeros([Nz, Nx]) # future time
c = 2000.0 #constant
# condition for not having dispersion (alias in space)
# d hale
fmax = c/(10.0*Ds)
fpeak = 0.5*fmax

# Koslov alias plane criteria (c must be the min value)
# fmax = c/(2.0*Ds)
# fpeak = 1.0*fmax
sourcewave = RickerSource(fpeak, 10.0, delay=1.0/fpeak)
sx = 0
sz = 0
V = np.zeros([Nz, Nx]) 
V[:][:] = c

# additional not needed 
Simulation  = np.zeros([numberiter, Nz, Nx])      


for i in range(0, numberiter):
    for k in range(Nz):
        for j in range(Nx):
            # u0k u1k*ujk*u3k u4k     
            # Boundary fixed 0 outside        
            u0k=u1k=u3k=u4k=0.0
            uj0=uj1=uj3=uj4=0.0
            ujk = Ucurrent[k][j]      

            if(j-2 > -1):
                u0k = Ucurrent[k][j-2]  
            if(j-1 > -1):
                u1k = Ucurrent[k][j-1]
            if(j+1 < Nx):
                u3k = Ucurrent[k][j+1]            
            if(j+2 < Nx):
                u4k = Ucurrent[k][j+2]
            if(k-2 > -1):
                uj0 = Ucurrent[k-2][j]
            if(k-1 > -1):
                uj1 = Ucurrent[k-1][j]
            if(k+1 < Nz):
                uj3 = Ucurrent[k+1][j]            
            if(k+2 < Nz):
                uj4 = Ucurrent[k+2][j]

            Ufuture[k][j] = 2*ujk-Uprevious[k][j]+(1.0/12.0)*(-u0k+16*u1k+16*u3k-u4k -uj0+16*uj1+16*uj3-uj4 -60*ujk)*(Dt*V[k][j]/Ds)**2

    # try adding the source to the laplacian

    # in the second loop you add the source to the laplacian
    # then finally you make the step forward

    # then finally  and then calculate the next step.
    # you can precalculate all weights as well.
    # add implementation using density could be nice
    # madagascar tip paul sava uses linear 
    # interpolation around the source and just 1 as dsm=1??


    # source position center grid
    # smooth region around the center of the grid +3-3
    # exact analytical solution circular
    if(i*Dt < 2.0/sourcewave.fc):
        dsm = 2
        for nz in range(max(sz-dsm,0),min(sz+dsm+1, Nz)):
            for nx in range(max(sx-dsm,0),min(sx+dsm+1, Nx)):
                dz = nz - sz 
                dx = nx - sx
                t = i*Dt
                r = np.sqrt(dz**2+dx**2)
                if(r == 0):
                    Ufuture[nz][nx] -= sourcewave(t)
                else:
                    Ufuture[nz][nx] -= sourcewave(t-r/c)/(2*np.pi*r)

    # make the update in the time stack
    Uprevious = Ucurrent
    Ucurrent = Ufuture
    Simulation[i] = Ucurrent

    sys.stderr.write("\r %d" %(i) )