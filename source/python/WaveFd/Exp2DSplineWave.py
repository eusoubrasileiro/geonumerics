#!/usr/bin/python
import sys
sys.path.append('../../python');

import numpy as np
import pylab as py
from Filters import WindowHann
from Wavelet import *
import time


def FourierDerivative(f):
    N = np.size(f)
    n = np.arange(0,N)
    # df discrete differential operator
    df = np.complex(0,1)*py.fftshift(n-N/2)
    dfdt = py.ifft( df*py.fft(f) )  
    return py.real(dfdt)
    



class FourierWaveField2D:
    """
    Explicit 2D wave equation
    3 order centered in time
    Fourier derivatives in space
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
        
        # calculate derivatives using fourier?
        # derivatives in x and y direction
        Dx = np.zeros([self.N, self.M]) 
        Dy = np.zeros([self.N, self.M])
        # xline by xline
        for i in range(self.N):        
            Dx[i] = FourierDerivative(self.Pn[i]) # first derivative
            Dx[i] = FourierDerivative(Dx[i]/self.Rho[i]) # second with Rho

        self.Pn = self.Pn.transpose()
        self.Rho = self.Rho.transpose()
        Dy = Dy.transpose()

        # yline by yline
        for i in range(self.M):        
            Dy[i] = FourierDerivative(self.Pn[i]) # first derivative
            Dy[i] = FourierDerivative(Dy[i]/self.Rho[i]) # second with Rho

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
    
    
    
