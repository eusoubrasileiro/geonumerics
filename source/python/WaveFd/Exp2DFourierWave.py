#!/usr/bin/python
import numpy as np
from BaseWave2DField import BaseWave2DField
import pylab as py

def FourierDerivative(f):
    N = np.size(f)
    n = np.arange(0,N)
    # df discrete differential operator
    df = np.complex(0,1)*py.fftshift(n-N/2)
    dfdt = py.ifft( df*py.fft(f) )  
    return py.real(dfdt)


class Exp2DFourierWave(BaseWave2DField):
    """
    Explicit 2D wave equation with 2 order centered in time
    and using spline for derivatives in space
        
    """

    def __init__(self,
                nx, 
                nz,
                ds,
                dt,
                velocity,
                rho,
                sx,
                sz,
                maxiter,
                nrec=1,
                wavelet=None):
        """
        Initialize a new wave equation field explicit centered differences time
        and spline for space.
        
        nx       : number of discretization in x
        nz       : number of discretization in z
        ds       : dx=dz=ds grid spacing 
        dt       : time step - e.g. seconds
        velocity : 2d velocity distribution
        rho      : 2d density distribution
        sx/sz    : source wavelet position in indexes (i, k)
        maxiter  : total iterations
        nrec     : recording interval 1 equals time step 
        wavelet  : source wavelet function applied at the position (sx, sz)
                   must have sample rate equal to dt
        """
        # creates the constant density part first
        super(Exp2DFourierWave, self).__init__(nx, nz, ds, dt, velocity, sx, sz, maxiter, nrec, wavelet)
        
        # centered differences in time, add previous time
        self.Uprevious = np.zeros([self.Nz, self.Nx])
        # Rho field for each (x, z) point
        # Rho doesnt vary with time
        self.Rho = np.zeros([self.Nz, self.Nx])
        # if a matrix of velocity is passed fills it        
        # setts up the velocity field
        if(type(rho) is np.ndarray):
            if(np.shape(rho) == (self.Nz, self.Nx) ):
                self.Rho = rho
        # if not put a constant velocity
        else:
            self.Rho[:][:] = rho


    def SolveNextTime(self):

        try:
            self.tstep += 1
        except :
            self.tstep = 0

        self.Source(self.tstep)
        
        # calculate derivatives using fourier
        # derivatives in x and y direction
        d2dx = np.zeros([self.Nz, self.Nx]) 
        d2dz = np.zeros([self.Nz, self.Nx])

        # xline by xline
        for i in range(self.Nz):        
            d2dx[i] = FourierDerivative(self.Ucurrent[i]) # first derivative
            d2dx[i] = FourierDerivative(d2dx[i]/self.Rho[i]) # second with Rho

        # to easy calculate on the z direction, transpose the matrixes
        self.Ucurrent = self.Ucurrent.transpose()
        self.Rho = self.Rho.transpose()
        d2dz = d2dz.transpose()
        # zline by zline
        for i in range(self.Nx):        
            d2dz[i] = FourierDerivative(self.Ucurrent[i]) # first derivative
            d2dz[i] = FourierDerivative(d2dz[i]/self.Rho[i]) # second with Rho
        # transpose back
        self.Ucurrent = self.Ucurrent.transpose()
        self.Rho = self.Rho.transpose()
        d2dz = d2dz.transpose()

        # derivatives combination + source
        LP = d2dx + d2dz      
        # simple centered differences in time 2nd order
        # matrix operations
        self.Ufuture = ((self.Vel*self.Dt)**2)*LP*self.Rho+2*self.Ucurrent-self.Uprevious     

        # update stack times
        self.Uprevious = self.Ucurrent
        self.Ucurrent = self.Ufuture        
        
        return self

    
    
