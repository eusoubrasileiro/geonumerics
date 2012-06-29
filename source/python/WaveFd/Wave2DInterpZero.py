#!/usr/bin/python
import sys
sys.path.append('../../python');

import numpy as np
import pylab as py
from Filters import SincLowPass
from Filters import WindowHann, WindowTukey
from Wavelet import *
import time
from numpy.polynomial import chebyshev as cheby
from scipy import interpolate as interp

"""
Tenth tentative using polynomial interpolations
to approximate the derivatives
now uses the 0 padding before calculating the derivatives.
"""

def FourierDerivative(f):
    """
    this derivatie just works for periodic 2*pi multiple series
    have to figure out how to make that work for any function
    """
    N = np.size(f)
    n = np.arange(0,N)
    # df discrete differential operator
    df = np.complex(0,1)*py.fftshift(n-N/2)
    dfdt = py.ifft( df*py.fft(f) )  
    return py.real(dfdt)
    


# amazing aproximation but very timing consuming
# and its also limited to the number of points that you use
def ChebyshevDerivative(x, y, m=4):
    """
    Fits a Chebyshev poly of degree (deg)
    in (y) and calculate the 1st derivative.
    
    (2*m+1) points are used to calculate the derivative
    in the borders previous/after points are used.
    """
    deg = 2*m # we always need 2*m+1 points for a degree 2*m
    n = np.size(x)
    if(n != np.size(y)
        or n < deg 
        or n < m*2+1):
        raise Exception(" size/degree/window inconsistence  ")

    # Chebyshev polys in the interval
    npoly = n-m*2
    Cheb_Coef = np.zeros([npoly, deg+1])          
    for i in range(npoly):
        Cheb_Coef[i] = cheby.chebfit(x[i:i+2*m+1], y[i:i+2*m+1], deg)

    # der coeficients    
    Cheb_Coefder = np.zeros([npoly, deg])          
    for i in range(npoly):
        Cheb_Coefder[i] = cheby.chebder(Cheb_Coef[i], 1)
   
    dydx = np.zeros(n)    
    dydx[0:m] = cheby.chebval(x[0:m], Cheb_Coefder[0]) #begin
    dydx[-m:] = cheby.chebval(x[-m:], Cheb_Coefder[npoly-1]) #end    
    for i in range(m, n-m): #middle
        dydx[i] = cheby.chebval(x[i], Cheb_Coefder[i-m])  

    return dydx

def SplineDerivative(x, y, k=3):
    """
    Fits a spline of degree (k)
    in (y) and calculate the 1st derivative.
    """
    n = np.size(x)
    if(n != np.size(y)
        or n < k):
        raise Exception(" size/degree inconsistence  ")

    Spline_Interp = interp.InterpolatedUnivariateSpline(x, y, k=k)
    dydx = np.zeros(np.size(x))
    # first derivative [1]
    for i in range(np.size(x)):
        dydx[i] = Spline_Interp.derivatives(x[i])[1]

    return dydx

def SplineDerivativeZero(Ds, y, k=3):
    """
    padd k zeros in the borders and calculate the derivative
    """
    # padd with zeros before and after size k
    padd = k+3
    y_ = np.append(np.zeros(padd), np.append(y, np.zeros(padd)))
    x = np.arange(0.0, Ds*np.size(y_), Ds)
    Spline_Interp = interp.InterpolatedUnivariateSpline(x, y_, k=k)
    dydx = np.zeros(np.size(y_))
    # first derivative [1]
    for i in range(np.size(x)):
        dydx[i] = Spline_Interp.derivatives(x[i])[1]

    return (dydx[padd:])[:-padd]

class WaveField2D:
    """
    Explicit 2D wave equation
    3 order centered in time
    Any derivatives in space
    """

    def __init__(self,
                 Ds=None,
                 Wavelet=None,
                 N=256,
                 M=256,
                 Dtr=0.04,
                 Si=50,
                 Sj=50,
                 Maxtime=0.5):
        """
        initialize a new wave equation field,
        for solving with finite differences method
        Ds = space increment in x (will be automatic modified if convergence can not be achieved based on Wavelet)
        Wavelet = source energy wavelet instance
        N number of discrete intervals in x  - ground dimension (e.g. meters)
        M
        Dtr time step for recording (e.g. seconds) 
        Si = energy source position
        Maxtime = simulation max time (seconds)
        Also use a better first approximation time to avoid bad wavelet formation
        """
        
        self.Ds = Ds
        self.N = N
        self.M = M
        self.Dtr = Dtr
        self.Si = Si
        self.Sj = Sj
        self.Maxtime = Maxtime
        
        if (isinstance(Wavelet, SourceWavelet) == False):
            raise Exception("Not a Wavelet type class")
        
        self._WvInst = Wavelet
        self.Wavelet = None
        
        # time step of solution, to be defined        
        self.Dt = None
        # amplitude or strain values, for each (x, z) point
        # centered in time
        self.Pn_0 = np.zeros([self.N, self.M])
        self.Pn = np.zeros([self.N, self.M])
        self.Pn_1 = np.zeros([self.N, self.M])
        # N is like line (i)        
        # M is like collum (j)
        # grid[N][M] = 0.0
        # velocity field for each (x) point
        # velocity doesn't vary with time
        self.Ve = None
        self.Rho = None

    def SetVe(self, Velocity):
        """
        sets the velocity field
        """
        # if a matrix of velocity is passed fills it
        if(type(Velocity) is np.ndarray):
            if(np.shape(Velocity) == [N, M] ):
                self.Ve = Velocity
        # if not put a constant velocity
        else:
            self.Ve = np.zeros([self.N, self.M]) + Velocity

    def SetRh(self, Density):
        """
        sets the Density field
        """
        # if a matrix of velocity is passed fills it
        if(type(Density) is np.ndarray):
            if(np.shape(Density) == [N, M] ):
                self.Rho = Velocity
        # if not put a constant density
        else:
            self.Rho = np.zeros([self.N, self.M]) + Density

    def Next(self):
        if(self.Ve == None 
           or self.Wavelet == None 
           or self.Dt == None
           or self.t == None):
            print "can't solve next steps"
            return
        
        # calculate derivatives using cheby
        # derivatives in x and y direction
        Dx = np.zeros([self.N, self.M]) 
        Dy = np.zeros([self.N, self.M])

        # chebyshev Good but very time consuming
        # number of points needed for the interpolation
        # to consider the wavelength Chebyshev
        #maxL = self.Ve.max()/self._WvInst.Fc
        #m = int(maxL/self.Ds)*2
        # ChebyshevDerivative(x, self.Pn[i], m)

        # set the borders as zero, tukey bichado and cos too
        # self.Pn[:] = WindowTukey(self.N)
        # self.Pn.transpose()[:] *= WindowTukey(self.M)    
        # self.Pn = self.Pn.transpose()

        # xline by xline
        for i in range(self.N):        
            Dx[i] = SplineDerivativeZero(self.Ds, self.Pn[i]) # first derivative
            Dx[i] = SplineDerivativeZero(self.Ds, Dx[i]/self.Rho[i]) # second with Rho
    
        self.Pn = self.Pn.transpose()
        self.Rho = self.Rho.transpose()
        Dy = Dy.transpose()

        # yline by yline
        for i in range(self.M):        
            Dy[i] = SplineDerivativeZero(self.Ds, self.Pn[i]) # first derivative
            Dy[i] = SplineDerivativeZero(self.Ds, Dy[i]/self.Rho[i]) # second with Rho
    
        self.Pn = self.Pn.transpose()
        self.Rho = self.Rho.transpose()
        Dy = Dy.transpose()

        # derivatives combination + source
        LP = Dx + Dy
        if(self.t < np.size(self.Wavelet)):
            LP[self.Si][self.Sj] -= self.Wavelet[self.t]
       
        # simple centered differences in time 2nd order
        self.Pn_1 = (self.Ve**2)*LP*self.Rho*self.Dt**2+2*self.Pn-self.Pn_0     

        # update stack times
        self.Pn_0 = self.Pn
        self.Pn = self.Pn_1        
        self.t = self.t + 1
        
        return self
    
    def _CurrentTime(self):
        if(self.Dt == None
           or self.t == None):
            print "can't solve next steps"
            
        return self.t*self.Dt
    
    def _GetDeltaSpace(self):
        """
        using Niquest principle for avoiding alias in space
        calculate Ds also using velocity = lambda * frequency 
        """
        if(self.Ve == None):
            raise Exception("Velocity field not set")
        
        Omega_fct = 0.2 # must be smaller than 1 to make the inequality true
        Vmin = self.Ve.min() # minimum  velocity
        Ds = Omega_fct*Vmin/(self._WvInst.Fc*2)
        
        # if required for convergence change it
        if(Ds < self.Ds or self.Ds == None):
            self.Ds = Ds
        
        return Ds 
        
       
    def _GetDeltaTime(self):
        """
        Calculate time step based on convergence criteria for 2D
        """
        J_fct = 0.1 # must be smaller than 1 to make the inequality true
        # remember time is 2 order so must be smaller for better convergence
        Ds = self._GetDeltaSpace()
        Vmax = self.Ve.max()
        self.Dt = 2*Ds*J_fct/(Vmax*np.pi)
        
        if(self.Dt > self.Dtr): # in a extreme case where recording step is smaller
            self.Dt = self.Dtr
            
        return self.Dt
    
    def _GetWavelet(self):
        """
        using the defined time step get the wavelet
        """
        if(self.Dt == None):
            raise Exception("Time step not defined")
            return
        
        self.Wavelet = self._WvInst.Samples(self.Dt)

    def _NumberInteractions(self):
        """
        Number of interaction based on Maxtime
        and time step
        """
        return int(self.Maxtime/self.Dt)
    
    def _IntervalInteractions(self):
        """
        number of interactions at every Dtr seconds
        """
        return int(self.Dtr/self.Dt) # snapshots interval at every Dtr seconds

    def TotalEstimatedTime(self):
        """
        Estimate time based on time for 100 interactions
        """
        if(self.Ve == None):
            raise Exception("Velocity field not set")
        
        self._GetDeltaTime()
        self._GetWavelet()
        
        initial = time.clock()
        self.t = 1

        for i in range(100):
            self.Next()

        final = time.clock()
        timeperstep = (final-initial)/100.0
        
        print "time per step (s)", timeperstep
        print "Estimated time (s)", self._NumberInteractions()*timeperstep
        print "Number of snapshots ", int(self._NumberInteractions()/self._IntervalInteractions())
        
        self.Pn_0[:][:] = 0.0 
        self.Pn[:][:] = 0.0
        self.Pn_1[:][:] = 0.0


    def Loop(self, Save=False, name='Exp1D'):
        """
        Loop through all time steps until (Maxtime)
        saving the matrix snapshots at every (Snapshots)
        """
        
        if(self.Ve == None):
            raise Exception("Velocity field not set")
        
        self._GetDeltaTime()
        self._GetWavelet()
        
        #raise        
        self.t = 1
        mt = self._NumberInteractions()
        nrc = self._IntervalInteractions() # snapshots interval at every Dtr seconds
        movie = np.zeros(((mt/nrc)+1, self.N)) # animation matrix at every nrc steps
        
        j=0 # counter for the animation
        initial = time.clock()
        
        for t in range(mt):
            self._Source()
            self.Next()
            if ( t%nrc == 0 ):
                movie[j] = self.Utime[1]
                j+=1
                
        final = time.clock()
         
        if(Save==True):
            np.save(name, movie)
        
        print "real total time (s) ", final-initial
        return movie
    
    def LoopReceiver(self, RcPos=10, Save=False, name='Exp1D'):
        """
        Loop through all time steps until (Maxtime)
        saving the values at the position RcPos 
        snapshots at every (Snapshots)
        """
        
        if(self.Ve == None):
            raise Exception("Velocity field not set")
        
        self._GetDeltaTime()
        self._GetWavelet()
        
        #raise        
        self.t = 1
        mt = self._NumberInteractions()
        nrc = self._IntervalInteractions() # snapshots interval at every Dtr seconds
        movie = np.zeros(((mt/nrc)+1)) # animation matrix at every nrc steps
        
        j=0 # counter for the animation
        initial = time.clock()
        
        for t in range(mt):
            self._Source()
            self.Next()
            if ( t%nrc == 0 ):
                movie[j] = self.Utime[1][RcPos]
                j+=1
                
        final = time.clock()
         
        if(Save==True):
            np.save(name, movie)
        
        print "real total time (s) ", final-initial
        return movie
    
    
    
