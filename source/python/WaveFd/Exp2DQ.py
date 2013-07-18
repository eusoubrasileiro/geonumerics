r"""

4th order explicit space test question

Convergence
 
.. math::
 
     \Delta t \leq \frac{2 \Delta s}{ V \sqrt{\sum_{a=-N}^{N} (|w_a^1| + |w_a^2|)}}


Where w_a are the centered differences weights

.. math::

    U_{jk}^{n+1}  =  \left( \frac{\Delta t  V_{jk}}{\Delta s} \right) ^2  \left[ \sum_{a=-N}^N  w_a \left( U_{j+a k}^n + U_{j k+a}^n \right) + S_{jk}^n {\Delta s}^2\right] + 2 U_{jk}^{n} - U_{jk}^{n-1}


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

Nz = Nx = 50
Dt = 0.003
Ds = 10.0
c = 2000.0 #constant
sx = 0
sz = 0
numberiter = 300
# condition for not having dispersion (alias in space)
# d hale
fmax = c/(10.0*Ds)
fpeak = 0.5*fmax
# Koslov alias plane criteria (c must be the min value)
# fmax = c/(2.0*Ds)
# fpeak = 1.0*fmax
sourcewave = RickerSource(fpeak, 1.0, delay=1.0/fpeak)
Uprevious = np.zeros([Nz, Nx]) # previous time
Ucurrent = np.zeros([Nz, Nx]) # current time
Ufuture = np.zeros([Nz, Nx]) # future time
Laplace = np.zeros([Nz, Nx]) # laplace operator or star convolution
Simulation  = np.zeros([numberiter, Nz, Nx])      
Vjk = np.zeros([Nz, Nx]) 

Vjk[:] = c
R2 = (Dt*Vjk/Ds)**2

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

            Laplace[k][j] = (1.0/12.0)*(-u0k+16*u1k+16*u3k-u4k -uj0+16*uj1+16*uj3-uj4 -60*ujk)

    # makes solution prettier!! 
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
                if(r == 0): # must be negative to propagate the exact wavelet.
                    Laplace[nz][nx] -= (Ds**2)*sourcewave(t)
                else:
                    Laplace[nz][nx] -= (Ds**2)*sourcewave(t-r/c)/(2*np.pi*r)

    Ufuture = 2*Ucurrent-Uprevious+Laplace*R2

    # make the update in the time stack
    Uprevious = Ucurrent
    Ucurrent = Ufuture
    Simulation[i] = Ucurrent

    sys.stderr.write("\r %d" %(i) )