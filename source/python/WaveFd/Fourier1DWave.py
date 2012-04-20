#!/usr/bin

import numpy as np
import time
import scipy.linalg as ln


class Wave1DField:
    """
    Implicit 1D wave equation
    2 order centered in space
    2 order backward in time Lax-Wendroff
    """

    def __init__(self,
                 Ds=None,
                 Wavelet=None,
                 N=200,
                 Dtr=0.01,
                 Si=100,
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
        self.Dt = 0.001
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
        self.Solved = False
        self.t = 1

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
        
    def R2(self, i):
        return (self.Ds/(self.Dt*self.Vel[i]))**2

        
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
        
        #if( t - 1 >= np.size(Wavelet)):
            # we should not force the solution after final time of wavelet
            # self.Utime[0][Si] = self.Utime[1][Si] = 0
            # return

        if(t - 1 < np.size(Wavelet)):
            self.Utime[0][Si] = Wavelet[t-1]
                
        if(t < np.size(Wavelet)):
            self.Utime[1][Si] = Wavelet[t]

        return
    
    def LinearSystem(self):
        """
        Assembly linear system
        """
        
        # assembly matrix of linear system
        # to solve u(t) based on u(t-1) and u(t-2)
        # the matrix includes all future values of u
        # in the entire grid, so size is the number of cells
        # start with zeros that is also the countour condition u(t)=0
        self.mUt = np.zeros([self.N, self.N])
        

        # assembly linear system, the linear system
        # ignores external part of the grid = locked boundary
        # ln go through all the cells in the grid Ut
        # each cell gives one equation (line)
        for Ln in range(0, self.N, 1): 
            # 1.0*u(x-1) + gama(x)*u(x) + 1.0*u(x+1) 
            # turn the indices to the one of original matrix
            i = Ln
            
            self.mUt[i][i] = -2.0-self.R2(i)

            if(i-1 >= 0): # u(x-1) inside grid in I
                self.mUt[Ln][Ln-1] = 1.0
            if(i+1 < self.N): # u(x+1) inside grid in I
                self.mUt[Ln][Ln+1] = 1.0
            
        return self.mUt

    def Independent(self):
        """
        Independent term
        """
        #independent term, where the previous times goes in
        self.vId = np.zeros([self.N])
        # fill the independent vector
        for Ln in range(0, self.N, 1): 
            # turn the indices to the one of original matrix
            i = Ln
            
            # -2 u(x,t-1) + u(x,t-2)
            self.vId[Ln] = -2*self.Utime[1][i]+self.Utime[0][i]
            # / (Ds/(Dt*Vel(x,z)))**2
            self.vId[Ln] *= self.R2(i)

        return self.vId


    def SolveSystem(self):
        """
        Find ... factorization of the matrix
        once found each time step is appying a recipe
        """
        self.LinearSystem()
        self.mUtfactor = ln.lu_factor(self.mUt)
        self.Solved = True

        return self.mUtfactor
        
        
    def Next(self):
        """
        Calculate the next time (factorization)
        and update the time stack grids
        """

        if(self.Solved == False):
            self.SolveSystem()
            self.t=1

        self._Source()
        # in time
        # As t is in [0, 1, 2] (2nd order)
        # time t in this case is Utime[2]

        v = self.Independent()

        result = ln.lu_solve(self.mUtfactor, v)
        # reshape the vector to became a matrix again
        self.Utime[2] = result

        # make the update in the time stack
        # before [t-2, t-1,  t]
        # after  [t-1, t,  t+1]
        # so t-2 receive t-1 and etc.

        u = self.Utime 
        u[0] = u[1]
        u[1] = u[2]
        
        self.t +=1

        return
    
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
        if(self.Vel == None):
            raise Exception("Velocity field not set")
        
        Omega_fct = 0.25 # must be smaller than 1 to make the inequality true
        Vmin = min(self.Vel) # minimum  velocity
        Ds = Omega_fct*Vmin/(self._WvInst.Fc*2)
        
        # if required for convergence change it
        if(Ds < self.Ds or self.Ds == None):
            self.Ds = Ds
        
        return Ds
    
    def _GetDeltaTime(self):
        """
        Calculate time step
        """
        
        if(self.Dt > self.Dtr or self.Dt == None): 
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
        if(self.Vel == None):
            raise Exception("Velocity field not set")
        
        self._GetDeltaTime()
        self._GetDeltaSpace()
        #self._GetWavelet()
        
        
        initial = time.clock()
        self.t = 1

        for i in range(100):
            self._Source()
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
        
        self._GetDeltaTime()
        self._GetDeltaSpace()
        #self._GetWavelet()
        
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
        
        if(self.Vel == None):
            raise Exception("Velocity field not set")
        
        self._GetDeltaSpace()()
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
    
    
