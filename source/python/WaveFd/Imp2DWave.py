#!/usr/bin
# import filters.py

#backend = 'gtk'
#import matplotlib
#matplotlib.use(backend)


import numpy as np
#from Filters import SincLowPass
# since the matrix is simetric Hermitian and positive definite
# we can use cholesky to 
import scipy.linalg as ln

"""
 Note: the wavelet must have the same sample rate 

def SincWavelet(N=101, Fc=40, dt=0.0005, plot=False):
    """
    Note: the wavelet must have the same sample rate 
    than the simulation!
    Create a wavelet (energy source) for the simulation.
"""
    N is number of samples for the wavelet
    fc is the F central frequency
    dt is the sample rate
    Sample the Box Sinc filter simetrically around the zero
    apply hanning window
    """
    print "total wavelet time : %.1f miliseconds" % (dt*N*1000)
    wavelet = SincLowPass(N, Fc, dt)
    #wavelet = wavelet*filters.WindowHann(N)
    #normalize [0,1]
    return wavelet/(np.max(wavelet)-np.min(wavelet))


def RickerWavelet(N=101, sg=0.5, dt=0.05):
    t = np.arange(-dt*(N-1)/2,(dt*(N-1)/2)+dt, dt)
    wv = 2/((3*sg)**0.5*np.pi**0.25)
    wv *= (1 - (t/sg)**2)*np.exp(-t**2/(2*sg**2))
    return wv


class WaveField:
    """
    Implicit wave equation (acoustic) , finite differences
    2 order centered in space
    2 order backward in time
    no convergence limitations??
    Change to Crank-Nicholson
    """

    def __init__(self,
                 Nx=100,
                 Nz=50,
                 Ds=0.5,
                 Dt=0.001,
                 Sx=10,
                 Sz=10,
                 Wavelet=None,
                 Fw=40,
                 Snapshots=1,
                 MaxIter=1000):
        """
        initialize a new wave equation field,
        for solving with finite diferences method
        Nx number of discretization in x  - ground dimension (e.g. meters)
        Nz number of discretization in z
        Ds = Dx = Dz grid spacing in x = grid spacing in z
        Dt time step (e.g. seconds) 
        (TODO: calculate criteria for convergence!!)
        (Sx, Sz) = energy source position
        Wavelet = wavelet position
        Snapshots = number of iterations between intervals
        TotalIter = number total of iterations
        TODO: dt has always to be much smaller than the desired
        time snapshots, and equal the wavelet sample rate
        use a variable for that after...
        """
        if(Wavelet == None):
            self.Wavelet=10*RickerWavelet(101,sg=Dt/10.0,dt=Dt) # 10 power

        self.Ds = Ds
        self.Dt = Dt
        self.Snapshots = Snapshots
        self.MaxIter = MaxIter
        self.Sx = Sx
        self.Sz = Sz
        
        # 2rd order on space, centered
        # due 2rd order finite diferences ... N+2 + order
        # countour definitions will be adressed when assembling the 
        # linear system
        self.Nx=Nx
        self.Nz=Nz
        # 2nd order time, backward
        # so we need plus order+2 grids
        self.Nt = 3
        
        # amplitude or strain values, for each (x, z) point
        # must follow the order nt, nz, nx
        self.Utime = np.zeros([self.Nt, self.Nz, self.Nx])
        # nx is like collums (i)
        # nz is like lines (k)
        # nt is like 4th dimension (t)
        # eg for firt time grid, 2nd line and 3rd colum
        # grid[0][1][2] = 0.0
        # grid[t][k][i] = 0.0
        # velocity field for each (x, z) point
        # velocity doesnt vary with time
        self.Vel = np.zeros([self.Nz, self.Nx])
        # linear system not solved yet
        self.Solved = False


    def SetVel(self, Velocity):
        """
        sets the velocy field
        """
        # if a matrix of velocity is passed fills it
        if(type(Velocity) is np.ndarray):
            if(np.shape(Velocity) == (self.Nz, self.Nx) ):
                self.Vel = Velocity
        # if not put a constant velocity
        else:
            self.Vel[:][:] = Velocity

    def alfa(self, k, i):
        return 0.5*(self.Vel[k][i]*self.Dt/self.Ds)**2


    def LinearSystem(self):
        """
        Assembly linear system
        """
        
        # assembly matrix of linear system
        # to solve u(t) based on u(t-1) and u(t-2)
        # the matrix includes all future values of u
        # in the entire grid, so size is the number of cells
        # start with zeros that is also the countour condition u(t)=0
        self.mUt = np.zeros([self.Nz*self.Nx, self.Nz*self.Nx])
        

        # assembly linear system, the linear system
        # ignores external part of the grid = locked boundary
        # ln go through all the cells in the grid Ut
        # each cell gives one equation (line)
        for Ln in range(0, self.Nz*self.Nx, 1): 
            # 1.0*u(x-1,z) + gama(x,z)*u(x,z) + 1.0*u(x+1,z) + 1.0*u(x,z-1) + 1.0*u(x,z+1) 
            # turn the indices to the one of original matrix
            i = Ln%self.Nx 
            k = Ln/self.Nx  

            self.mUt[Ln][Ln] = 4*self.alfa(k, i)+1
            #is this right?
            if(i-1 >= 0): # u(x-1,z) inside grid in I
                self.mUt[Ln][Ln-1] = -self.alfa(k, i)
            if(i+1 < self.Nx): # u(x+1,z) inside grid in I
                self.mUt[Ln][Ln+1] = -self.alfa(k, i)
            if(k-1 >= 0): #u(x,z-1)
                self.mUt[Ln][Ln-self.Nx]= -self.alfa(k, i)
            if(k+1 < self.Nz): #u(x,z+1)
                self.mUt[Ln][Ln+self.Nx]= -self.alfa(k, i)

        return self.mUt

    def Independent(self):
        """
        Independent term
        """
        #independent term, where the previous times goes in
        self.vId = np.zeros([self.Nz*self.Nx])
        u = self.Utime
        # fill the independent vector
        for Ln in range(0, self.Nz*self.Nx, 1): 
            # turn the indices to the one of original matrix
            i = Ln%self.Nx 
            k = Ln/self.Nx         
            # boundary locked    
            u0 = u1 = u2 = u3 = 0.0
            
            if(i-1 >= 0): # u(x-1,z) inside grid in I
                u0 = u[1][k][i-1]
            if(i+1 < self.Nx): # u(x+1,z) inside grid in I
                u1 = u[1][k][i+1]
            if(k-1 >= 0): #u(x,z-1)
                u2 = u[1][k-1][i]
            if(k+1 < self.Nz): #u(x,z+1)
                u3 = u[1][k+1][i]
            
            
            self.vId[Ln] = self.alfa(k, i)*(u0+u1+u2+u3)
            self.vId[Ln] += (2-4*self.alfa(k, i))*u[1][k][i] - u[0][k][i]

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
        
    # problems may arise if you dont set
    # the system initial boundary condition
    # before assemblying and solving
    # the system. Also should be good
    # compare the difference matrices factors
    # in different stages of the solution
    # to see if the lu solution is really convergin...
    def SolveNextTime(self):
        """
        Calculate the next time (factorization)
        and update the time stack grids
        """

        if(self.Solved == False):
            self.SolveSystem()
            self.iter=1
            self.SourceBoundaryCondition(self.iter)            
            self.SolveSystem()
        
        self.SourceBoundaryCondition(self.iter)
        # in time
        # As t is in [0, 1, 2] (2nd order)
        # time t in this case is Utime[2]

        v = self.Independent()

        result = ln.lu_solve(self.mUtfactor, v)
        # reshape the vector to became a matrix again
        self.Utime[2] = np.reshape(result, (self.Nz, self.Nx))

        # make the update in the time stack
        # before [t-2, t-1,  t]
        # after  [t-1, t,  t+1]
        # so t-2 receive t-1 and etc.

        u = self.Utime 
        u[0] = u[1]
        u[1] = u[2]
        
        self.iter +=1

        return
    

    def SourceBoundaryCondition(self, it):
        """
        ( wavelet ) Set the boundary condition at the pertubation source position.
        ( sx, sz ) source position
        ( it ) At the given time step. Set t and t-1. (2 order finite diferences)
        if the iteration time is greater than the wavelet time
        sets the source position as 0
        TODO:
        verify if its really working
        """
        Wavelet=self.Wavelet
        Sx=self.Sx
        Sz=self.Sz
        
        if( it - 1 >= np.size(Wavelet)):
            self.Utime[0][Sz][Sx] = self.Utime[1][Sz][Sx] = 0
            return

        if(it - 1 < np.size(Wavelet)):
            self.Utime[0][Sz][Sx] = Wavelet[it-1]
                
        if(it < np.size(Wavelet)):
            self.Utime[1][Sz][Sx] = Wavelet[it]

        return

    def Loop(self, Save=True):
        """
        Loop through all time steps until (MaxIter)
        saving the matrix snapshots at every (Snapshots)
        """
        MaxIter=self.MaxIter 
        Snapshots=self.Snapshots         
        
        # for little problems with the wavelet put initialize as 1
        for i in range(1, MaxIter, 1):
            self.SourceBoundaryCondition(i)
            self.SolveNextTime()
            if(i%Snapshots==0): # every n'th Snapshots
                np.save("IfE"+str(i), self.Utime[1])

        return
