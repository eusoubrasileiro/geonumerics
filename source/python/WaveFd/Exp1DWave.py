#!/usr/bin

import numpy as np

def LinearSin(Fc=40.0, dt=None, plot=False):
    """
    Linear decreasing one period sin(2pi*f) 
    """
    # frequency of niquest limit
    if ( dt > 1/(2.0*Fc) or dt == None):
        dt = 1/(2.0*Fc) 
    
    t = np.arange(0, 1/Fc, dt)
    wavelet = np.sin(2*np.pi*Fc*t)*(-Fc*t+1)
    
    print "total wavelet time : %.1f miliseconds" % (dt*np.size(wavelet)*1000)
    wavelet = wavelet/(np.max(wavelet)-np.min(wavelet))
    
    # invert the function so, it starts with a small perturbation
    return wavelet[::-1]

class Wave1DField:
    """
    Explicit 1D wave equation
    2 order centered in space
    4 order backward in time
    """

    def __init__(self,
                 N=100,
                 Ds=0.5,
                 Dtr=0.04,
                 Si=50,
                 Maxtime=2.0):
        """
        initialize a new wave equation field,
        for solving with finite diferences method
        Nx number of discretization in x  - ground dimension (e.g. meters)
        Ds = Dx grid spacing in x 
        Dtr time step for recording (e.g. seconds) 
        Si = energy source position
        Wavelet = source energy wavelet
        Maxtime = simulation max time (seconds)
        TODO: dt has always to be much smaller than the desired
        time snapshots, and equal the wavelet sample rate
        use a variable for that after...
        """
        self.N = N
        self.Ds = Ds
        self.Dtr = Dtr
        self.Si = Si
        self.Maxtime = Maxtime
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
        self.Vel = np.zeros([self.N])

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
            self.Vel[:][:] = Velocity

    def _Source(self):
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
        
        if( t - 1 >= np.size(Wavelet)):
            self.Utime[0][Si] = self.Utime[1][Si] = 0
            return

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
        
        self._Source()
        
        for j in range(self.N):
            
            # sequence
            # 0, 1, 2 => j-1, j, j+1
            # 0, 1, 2 => n-1, n, n+1
            ej1n0 = self.Utime[0][j] # n-1
            
            ej0n1=0.0
            ej1n1 = self.Utime[1][j] # n
            ej2n1=0.0
            
            ej1n2 = 0.0 # n+1
            
            if j-1 > 0: 
                ej0n1 = self.Utime[1][j-1]
            if j+1 < self.N: 
                ej2n1 = self.Utime[1][j+1]
                
            ej1n2 = (ej0n1-2*ej1n1+ej2n1)*(self.Dt*self.Vel[j]/self.Ds)**2
            ej1n2 += 2*ej1n1-ej1n0
            
            
            # Lax-Wendroff 4 order time correction Hjn
            # space derivatives to solve time derivatives 
            # sequence
            # 0, 1, 2, 3, 4 => j-2, j-1, j, j+1, j+2
            # 0, 1, 2 => n-1, n, n+1
            # there is no propagation outside boundaries
            ej0n = 0.0
            ej1n = 0.0
            ej2n = self.Utime[1][j]
            ej3n = 0.0
            ej4n = 0.0
            # constant velocity outside boundaries
            vj1n = self.Vel[j]
            vj2n = self.Vel[j]
            vj3n = self.Vel[j]
            
            if j-1 > 0: 
                ej1n = self.Utime[1][j-1]
                vj1n = self.Vel[j-1]
            if j+1 < self.N: 
                ej3n = self.Utime[1][j+1]
                vj3n = self.Vel[j+1]
            if j-2 > 0:
                ej0n = self.Utime[1][j-2]
            if j+2 < self.N:
                ej4n = self.Utime[1][j+2]
            
            Hjn = 2*vj2n*(vj1n-2*vj2n+vj3n)*(ej1n-2*ej2n+ej3n)/(self.Ds**4)
            Hjn += (vj2n**2)*(ej0n-4*ej1n+6*ej2n-4*ej3n+ej4n)/(self.Ds**4)
            Hjn *= (self.Dt**2)/12
            self.Utime[2][j] =  ej1n2 - (Hjn*self.Dt**2)
        
        # update stack times
        self.Utime[0] = self.Utime[1]
        self.Utime[1] = self.Utime[2]
        
        self.t = self.t + 1
        
        return self

    def Loop(self, Save=True):
        """
        Loop through all time steps until (Maxtime)
        saving the matrix snapshots at every (Snapshots)
        """
        
        if(self.Vel == None 
           or self.Wavelet == None 
           or self.Dt == None):
            print "can't solve next steps"
            return
        #raise        
        self.t = 1
        for t in range(int(self.Maxtime/self.Dt)):
            self._Source()
            self.Next()
            if ( t*self.Dt%self.Dtr == 0 ): # every Dtr seconds
                np.save("IfE"+str(t), self.Utime[1])

        return
    
    