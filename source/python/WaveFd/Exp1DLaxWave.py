r"""
Need review, to base class usage and new wavelet classes
include also convergence criteria
"""

import numpy as np
import time        
from Wavelet import SourceWavelet

class LaxWand1DWave:
    """
    Explicit 1D wave equation
    4 order centered in space
    4 order centered in time Lax-Wendroff
    """

    def __init__(self,
                 Ds=None,
                 Wavelet=None,
                 N=300,
                 Dtr=0.005,
                 Si=150,
                 Maxtime=0.5):
        """
        initialize a new wave equation field,
        for solving with finite differences method
        Ds = space increment in x (will be automatic modified if convergence can not be achieved based on Wavelet)
        Wavelet = source energy wavelet instance
        Nx number of discrete intervals in x  - ground dimension (e.g. meters)
        Dtr time step for recording (e.g. seconds) 
        Si = energy source position
        Maxtime = simulation max time (seconds)
        Also use a better first approximation time to avoid bad wavelet formation
        """
        
        self.Ds = Ds
        self.N = N
        self.Dtr = Dtr
        self.Si = Si
        self.Maxtime = Maxtime
        
        if (isinstance(Wavelet, SourceWavelet) == False):
            raise Exception("Not a Wavelet type class")
        
        self._WvInst = Wavelet
        self.Wavelet = None
        
        # time step of solution, to be defined        
        self.Dt = None
        # 2nd order time, backward
        # so we need plus order+2 grids
        self.Nt = 3
        # amplitude or strain values, for each (x) point
        # must follow the order nt, n
        self.Utime = np.zeros([self.Nt, self.N])
        # N is like collum (i)        
        # Nt is like 2nd dimension (t)
        # eg for first time line grid
        # grid[0][1] = 0.0
        # grid[t][i] = 0.0
        # velocity field for each (x) point
        # velocity doesn't vary with time
        self.Vel = None

    def SetVel(self, Velocity):
        """
        sets the velocy field
        """
        # if a matrix of velocity is passed fills it
        if(type(Velocity) is np.ndarray):
            if(np.shape(Velocity) == (self.N) ):
                self.Vel = Velocity
        # if not put a constant velocity
        else:
            self.Vel = np.zeros(self.N) + Velocity
            

        
    def Source(self):
        """
        ( Wavelet ) Set the boundary condition at the pertubation source position.
        ( sx ) source position
        ( t ) At the given time step. Set t and t-1. 
        if the iteration time is greater than the Wavelet time
        sets the source position as 0
        """
        t = self.t
        
        Wavelet=self.Wavelet
        Si=self.Si        
        
        #if( t - 1 >= np.size(Wavelet)):
            # we should not force the solution after final time of wavelet
            # self.Utime[0][Si] = self.Utime[1][Si] = 0
            # return

        if(t - 1 < np.size(Wavelet)):
            self.Utime[0][Si] = Wavelet[t-1]
                
        if(t < np.size(Wavelet)):
            self.Utime[1][Si] = Wavelet[t]

        return
    
    def Next(self):
        if(self.Vel == None 
           or self.Wavelet == None 
           or self.Dt == None
           or self.t == None):
            print "can't solve next steps"
            return
        
        self.Source()
        
        # for each spatial position
        for j in range(self.N):
            
            # sequence 4 order space, 4 order time
            # 0, 1, 2, 3, 4, 5, 6 => j-3, j-2, j-1, j, j+1, j+2, j+3
            # 0, 1, 2 => n-1, n, n+1
            uj3n0 = self.Utime[0][j] # n-1
            uj3n=self.Utime[1][j] # n
            uj3n2 = 0.0 # n+1			
            # just for convention n means n1 = time n
            # BOUNDARY CONDITION TO BE LOOKED AT!
            # no propagtion outside boundaries?			
            uj0n=0.0
            uj1n=0.0
            uj2n=0.0            
            uj4n=0.0
            uj5n=0.0
            uj6n=0.0      
            # Lax-Wendroff 4 order time correction LWjn
            # space derivatives to solve time derivatives 
            # sequence, space derivatives 4 order
            # 0, 1, 2, 3, 4, 5, 6 => j-3, j-2, j-1, j, j+1, j+2, j+3
            # 0, 1, 2 => n-1, n, n+1
            # 1, 2, 3, 4, 5, => v-2, v-1, v, v+1, v+2
            # BOUNDARY CONDITION TO BE LOOKED AT!
            # there is no propagation outside boundaries??
            # constant velocity outside boundaries???
            # whyyyyy??? boundary condition weird???
            vj1 = self.Vel[j]
            vj2 = self.Vel[j]
            vj3 = self.Vel[j]
            vj4 = self.Vel[j]
            vj5 = self.Vel[j]
            # velocity is just need until 4th order, that's why just 5 indexes
            LWjn=0.0
            #######################################
            
            if j-3 > 0:
                uj0n = self.Utime[1][j-3]
            if j-2 > 0: 
                uj1n = self.Utime[1][j-2]
                vj1 = self.Vel[j-2]
            if j-1 > 0: 
                uj2n = self.Utime[1][j-1]
                vj2 = self.Vel[j-1]
            if j+1 < self.N: 
                uj4n = self.Utime[1][j+1]
                vj4 = self.Vel[j+1]
            if j+2 < self.N: 
                uj5n = self.Utime[1][j+2]
                vj5 = self.Vel[j+2]
            if j+3 < self.N:
                uj6n = self.Utime[1][j+3]
            
            # simple forth order space, 2 order time
            uj3n2 = (-uj1n+16*uj2n-30*uj3n+16*uj4n-uj5n)/12
            uj3n2 *= (self.Dt*vj3)**2/(self.Ds**2)
            uj3n2 += 2*uj3n-uj3n0
            
            # Lax-Wendroff 4 order time correction LWjn
            LWjn = ((vj1-8*vj2+8*vj4-vj5)**2)/144
            LWjn += vj3*(-vj1+16*vj2-30*vj3+16*vj4-vj5)/12
            LWjn *= (-uj1n+16*uj2n-30*uj3n+16*uj4n-uj5n)/6
            LWjn += vj3*(vj1-8*vj2+8*vj4-vj5)*(uj0n-8*uj1n+13*uj2n-13*uj4n+8*uj5n-uj6n)/48
            LWjn += (vj3**2)*(-uj0n+12*uj1n-39*uj2n+56*uj3n-39*uj4n+12*uj5n-uj6n)/6
            LWjn *= ((self.Dt*vj3)**2)/(12*self.Ds**4)
            self.Utime[2][j] =  uj3n2 + LWjn*(self.Dt**2)
        
        # update stack times
        self.Utime[0] = self.Utime[1]
        self.Utime[1] = self.Utime[2]
        
        self.t = self.t + 1
        
        return self
    
    def CurrentTime(self):
        if(self.Dt == None
           or self.t == None):
            print "can't solve next steps"
            
        return self.t*self.Dt
    
    def GetDeltaSpace(self):
        """
        using Niquest principle for avoiding alias in space
        calculate Ds also using velocity = lambda * frequency 
        """
        if(self.Vel == None):
            raise Exception("Velocity field not set")
        
        Omega_fct = 0.5 # must be smaller than 1 to make the inequality true
        Vmin = min(self.Vel) # minimum  velocity
        Ds = Omega_fct*Vmin/(self._WvInst.Fc*2)
        
        # if required for convergence change it
        if(Ds < self.Ds or self.Ds == None):
            self.Ds = Ds
        
        return Ds 
        
    def CharacteristicR(self):
        """
        Convergence criteria for 1D
        Lax-Wendroff 4order time and space
        Look at Jing-Bo Chen Geophysics 
        R expression
        """
        a=64.0 # sum of modulus second spatial derivatives weights
        b=160.0 # sum of modulus forth spatial derivatives weights
        return 2*np.sqrt(6)/np.sqrt(3*a+np.sqrt(9*a**2+12*b))
        
    def GetDeltaTime(self):
        """
        Calculate time step based on convergence criteria for 1D
        Based on:
        1) Dan D. Kosloff and Edip Baysal (Forward modeling by a Fourier mehtod) - APPENDIX
        2) Panorama Technologies Chapter 2 (Finite Differences page 56) 
        """
        J_fct = 0.5 # must be smaller than 1 to make the inequality true
        Ds = self.GetDeltaSpace()
        Vmax = max(self.Vel)
        self.Dt = 2*Ds*J_fct/(Vmax*2*3.141592654)
        
        if(self.Dt > self.Dtr): # in a extreme case where recording step is smaller
            self.Dt = self.Dtr
            
        return self.Dt
    
    def GetWavelet(self):
        """
        using the defined time step get the wavelet
        """
        if(self.Dt == None):
            raise Exception("Time step not defined")
            return
        
        self.Wavelet = self._WvInst.Samples(self.Dt)

    def TotalEstimatedTime(self):
        """
        Estimate time based on time for 100 interactions
        """
        if(self.Vel == None):
            raise Exception("Velocity field not set")
        
        if(self.Ds == None):
            self.GetDeltaSpace()
        
        if(self.Dt == None):
            self.GetDeltaTime()
        
        initial = time.clock()
        self.t = 1

        for i in range(100):
            self.Source()
            self.Next()

        final = time.clock()
        timeperstep = (final-initial)/100.0
        
        print "time per step (s)", timeperstep
        print "Estimated time (s)", self._NumberInteractions()*timeperstep
        print "Number of snapshots ", int(self._NumberInteractions()/self._IntervalInteractions())
        
        self.Utime[:][:] = 0

    def Loop(self, Save=False, name='Exp1D'):
        """
        Loop through all time steps until (Maxtime)
        saving the matrix snapshots at every (Snapshots)
        """
        
        if(self.Vel == None):
            raise Exception("Velocity field not set")
        
        if(self.Ds == None):
            self.GetDeltaSpace()
        
        if(self.Dt == None):
            self.GetDeltaTime()
        
        self.GetWavelet()
        
        #raise        
        self.t = 1
        mt = int(self.Maxtime/self.Dt)
        nrc = int(self.Dtr/self.Dt) # snapshots interval at every Dtr seconds
        movie = np.zeros(((mt/nrc)+1, self.N)) # animation matrix at every nrc steps
        
        j=0 # counter for the animation
        initial = time.clock()
        
        for t in range(mt):
            self.Source()
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
        
        if(self.Vel == None):
            raise Exception("Velocity field not set")
        
        if(self.Ds == None):
            self.GetDeltaSpace()
        
        if(self.Dt == None):
            self.GetDeltaTime()
        
        self.GetWavelet()
        
        #raise        
        self.t = 1
        mt = int(self.Maxtime/self.Dt)
        nrc = int(self.Dtr/self.Dt) # snapshots interval at every Dtr seconds
        movie = np.zeros(((mt/nrc)+1)) # animation matrix at every nrc steps
        
        j=0 # counter for the animation
        initial = time.clock()
        
        for t in range(mt):
            self.Source()
            self.Next()
            if ( t%nrc == 0 ):
                movie[j] = self.Utime[1][RcPos]
                j+=1
                
        final = time.clock()
         
        if(Save==True):
            np.save(name, movie)
        
        print "real total time (s) ", final-initial
        return movie
    
    
