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

def _ricker(fc, dt, t, delay=True):
    if(delay):
        tdelay = 1.0/fc # to create a causal wavelet starting at 0
        t = t-tdelay
    
    ricker = (1-2*(np.pi*fc*t)**2)*np.exp(-(np.pi*fc*t)**2) 
    return  ricker

def Ricker(fc, dt, n=None):
    r"""
    Ricker Wavelet:
    :math:`A = (1-2 \pi^2 f^2 t^2) e^{-\pi^2 f^2 t^2}`

    * fc        : maximum desired frequency
    * dt        : sample rate
    """
    # n is specified based on closer to zero 
    # 0.0 = (1-2 pi^2 fc^2 t^2) e^(-pi^2 fc^2 t^2)
    # that gives {t = -/+ 1/(sqrt(2) pi fc)}
    
    if(n==None):
        n=int(1/float(fc*dt))

    t = np.arange(-dt*(n-1)/2,(dt*(n-1)/2)+dt, dt)
    
    return  _ricker(fc, dt, t)

Nz = Nx = 80
Dt = 0.002
Ds = 10
numberiter = 2000



#Source = Ricker(40, 0.001, 80)
Uprevious = np.zeros([Nz, Nx]) # previous time
Ucurrent = np.zeros([Nz, Nx]) # current time
Ufuture = np.zeros([Nz, Nx]) # future time
c = 2000.0 #constant
fmax = c/(10.0*Ds)
fpeak = 0.5*fmax
sx = 40
sz = 40
V = np.zeros([Nz, Nx]) 
V[:][:] = c

# additional not needed 
Simulation  = np.zeros([numberiter, Nz, Nx])      


for i in range(0, numberiter):
    for k in range(Nz):
        for j in range(Nx):
            # u0k u1k*ujk*u3k u4k     
            # Boundary fixed 0 outside        
            u0k=u1k=0.0
            uj0=uj1=0.0
            ujk = Ucurrent[k][j]      

            if(j-1 > -1):
                u0k = Ucurrent[k][j-1]
            if(j+1 < Nx):
                u1k = Ucurrent[k][j+1]            
            if(k-1 > -1):
                uj0 = Ucurrent[k-1][j]
            if(k+1 < Nz):
                uj1 = Ucurrent[k+1][j]            

            Ufuture[k][j] = 2*ujk-Uprevious[k][j]+(uj0+uj1+u0k+u1k-4*ujk)*(Dt*V[k][j]/Ds)**2

    # source position center grid
    # smooth region around the center of the grid +3-3
    # exact solution
    dsm = 3
    for nz in range(sz-dsm,sz+dsm,1):
        for nx in range(sx-dsm,sx+dsm,1):
            dz = nz - sz 
            dx = nx - sx
            t = i*Dt
            r = np.sqrt(dz**2+dx**2)
            if(r == 0):
                Ufuture[nz][nx] -= _ricker(fpeak, Dt, t-r/c, delay=True)
            else:
                Ufuture[nz][nx] -= _ricker(fpeak, Dt, t-r/c, delay=True)/(2*np.pi*r)

    # make the update in the time stack
    Uprevious = Ucurrent
    Ucurrent = Ufuture
    Simulation[i] = Ucurrent

    sys.stdout.write("\r %d" %(i) )