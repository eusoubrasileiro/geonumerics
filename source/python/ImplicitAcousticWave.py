#!/usr/bin
# import filters.py

#backend = 'gtk'
#import matplotlib
#matplotlib.use(backend)

# TODO
# since the velocity field
# space interval and increment in time
# don't change in time, the main matrix
# (linear system) never changes
# what changes is the independent term
# could use for example LU decomposition
# to increase performeance
# instead of solving in every step
# all the system again

import pylab as py
import numpy as np
import time
import os
import filters


def SincWavelet(N=101, Fc=40, dt=0.0005, plot=False):
    """
    Note: the wavelet must have the same sample rate 
    than the simulation!
    Create a wavelet (energy source) for the simulation.
    N is number of samples for the wavelet
    fc is the F central frequency
    dt is the sample rate
    Sample the Box Sinc filter simetrically around the zero
    apply hanning window
    """
    print "total wavelet time : %.1f miliseconds" % (dt*N*1000)
    wavelet = filters.SincLowPass(N, Fc, dt)
    #wavelet = wavelet*filters.WindowHann(N)
    #normalize [0,1]
    return wavelet/(np.max(wavelet)-np.min(wavelet))



    
__doc__ = """ Implicit wave equation (acoustic) , finite differences
              2 order centered in space
              2 order backward in time
              no convergence limitations??
          """

class WaveField:
    def __init__(self, Nx=100, Nz=50, Ds=0.5, Dt=0.001, 
                 Sx=2, Sz=2, Wavelet=SincWavelet(Fc=40,dt=0.001), 
                 Snapshots=1, MaxIter=1000):
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
        self.Wavelet=10*Wavelet # 10 power
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
        self.SetVel(2000)

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

    def gama(self, k, i):
        return -(4.0+(self.Ds/(self.Dt*self.Vel[k][i]))**2)


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
        # ignores external part of the grid = free boundary
        # ln go through all the cells in the grid Ut
        # each cell gives one equation (line)
        for Ln in range(0, self.Nz*self.Nx, 1): 
            # 1.0*u(x-1,z) + gama(x,z)*u(x,z) + 1.0*u(x+1,z) + 1.0*u(x,z-1) + 1.0*u(x,z+1) 
            # turn the indices to the one of original matrix
            i = Ln%self.Nx 
            k = Ln/self.Nx  

            self.mUt[Ln][Ln] = self.gama(k, i)

            if(i-1 > 0): # u(x-1,z) inside grid in I
                self.mUt[Ln][Ln-1] = 1.0
            if(i+1 < self.Nx): # u(x+1,z) inside grid in I
                self.mUt[Ln][Ln+1] = 1.0
            if(k-1 > 0): #u(x,z-1)
                self.mUt[Ln][Ln-self.Nx]= 1.0
            if(k+1 < self.Nz): #u(x,z+1)
                self.mUt[Ln][Ln+self.Nx]= 1.0

        return self.mUt

    def Independent(self):
        """
        Independent term
        """
         #independent term, where the previous times goes in
        self.vId = np.zeros([self.Nz*self.Nx])
        # fill the independent vector
        for Ln in range(0, self.Nz*self.Nx, 1): 
            # turn the indices to the one of original matrix
            i = Ln%self.Nx 
            k = Ln/self.Nx  
            
            # -2 u(x,z,t-1) + u(x,z,t-2)
            self.vId[Ln] = -2*self.Utime[1][k][i]+self.Utime[0][k][i]
            # / (Ds/(Dt*Vel(x,z)))**2
            self.vId[Ln] /= (self.Ds/(self.Dt*self.Vel[k][i]))**2

        return self.vId
    

    def SourceBoundaryCondition(self, it):
        """
        ( wavelet ) Set the boundary condition at the pertubation source position.
        ( sx, sz ) source position
        ( it ) At the given time step. Set t and t-1. (2 order finite diferences)
        if the iteration time is greater than the wavelet time
        sets the source position as 0
        """
        Wavelet=self.Wavelet
        Sx=self.Sx
        Sz=self.Sz
        
        if(it - 1 >= 0):
            if(it - 1 > np.size(Wavelet) ):
                self.Utime[0][Sz][Sx] = 0
                return
            else:
                self.Utime[0][Sz][Sx] = Wavelet[it-1]
                
        if(it > np.size(Wavelet)):
            self.Utime[1][Sz][Sx] = 0
        else:
            self.Utime[1][Sz][Sx] = Wavelet[it]
        
        return            
    
    
    def Loop(self, MakeGif=True):
        """
        Loop through all time steps until (MaxIter)
        taking picture snapshots at every (Snapshots)
        MakeGif : create a Gif animation with all the snapshots
        """
        MaxIter=self.MaxIter 
        Snapshots=self.Snapshots 
        
        py.ion()
        img = py.imshow(self.Utime[1])
        py.show()
         
        # for little problems with the wavelet put initialize as 1
        for i in range(1, MaxIter, 1):
            self.SourceBoundaryCondition(i)
            self.SolveNextTime()
            if(i%Snapshots==0): # every n'th Snapshots
                img.set_data(self.Utime[1])
                py.draw()
                time.sleep(0.1)
                if(MakeGif):
                    py.savefig("IfE"+str(i))
                
        if(MakeGif):
            #create gif animation
            os.system("convert -delay 50 -dispose None IfE*.png -loop 0 InfEn.gif")
            os.system("rm IfE*.png")
            
        return
        
    def SolveNextTime(self):
        """
        Calculate the next time
        and update the time stack grids
        """

        # in time
        # As t is in [0, 1, 2] (2nd order)
        # time t in this case is Utime[2]
        
        m = self.LinearSystem()
        v = self.Independent()
        result = np.linalg.solve(m, v)
        # reshape the vector to became a matrix again
        self.Utime[2] = np.reshape(result, (self.Nz, self.Nx))

        # make the update in the time stack
        # before [t-2, t-1,  t]
        # after  [t-1, t,  t+1]
        # so t-2 receive t-1 and etc.

        u = self.Utime 
        u[0] = u[1]
        u[1] = u[2]
        
      
        return


