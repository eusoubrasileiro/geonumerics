#!/usr/bin

import numpy as np
import scipy.linalg as ln
from Wavelet import Triangle
import gc, sys

class Imp2DLuWave:
    """
    Implicit wave equation (acoustic) , finite differences
    2 order centered in space
    2 order backward in time
    using Lu decomposition 
    """        
    def __init__(self,
                 nx,
                 nz,
                 ds,
                 dt,
                 velocity,
                 sx,
                 sz,
                 maxiter,
                 nrec=1,
                 wavelet=None,
                 ):
        """
        Initialize a new wave equation field implict centered differences second order
        
        nx       : number of discretization in x
        nz       : number of discretization in z
        ds       : dx=dz=ds grid spacing 
        dt       : time step - e.g. seconds
        velocity : 2d velocity distribution
        sx/sz    : source wavelet position in indexes (i, k)
        maxiter  : total iterations
        nrec     : recording interval 1 equals time step 
        wavelet  : source wavelet function applied at the position (sx, sz)
                   must have sample rate equal to dt
        """
        self.Ds = ds
        self.Dt = dt
        self.Nrec = nrec
        self.MaxIter = maxiter
        self.Si = sx
        self.Sk = sz
        self.Wavelet = wavelet
        
        # 2rd order on space, centered
        # due 2rd order finite diferences ... N+2 + order
        # countour definitions will be adressed when assembling the 
        # linear system
        self.Nx=nx
        self.Nz=nz
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
        # if a matrix of velocity is passed fills it
        
        # setts up the velocity field
        if(type(velocity) is np.ndarray):
            if(np.shape(velocity) == (self.Nz, self.Nx) ):
                self.Vel = velocity
        # if not put a constant velocity
        else:
            self.Vel[:][:] = velocity
    
        if(wavelet == None):         
            # using the principle of planar waves 
            # to set the source wavelet frequency               
            minvelocity = self.Vel[0][0]
            for array in self.Vel:
                vmin = array.min()
                if(minvelocity > vmin):
                    minvelocity = vmin                
            #default wavelet triangular
            self.Fw = minvelocity/(2*self.Ds);    
            self.Wavelet=10*Triangle(Fc=self.Fw,Dt=self.Dt)


    def Alpha(self, k, i):
        """
        repeats a lot in the code and in the matrix so
        better a def or it, if V varies
        """
        return 0.5*(self.Vel[k][i]*self.Dt/self.Ds)**2


    def LinearSystem(self):
        """
        Assembly linear system
        Depends on Velocity field and Alpha
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

            self.mUt[Ln][Ln] = 4*self.Alpha(k, i)+1
            #is this right?
            if(i-1 >= 0): # u(x-1,z) inside grid in I
                self.mUt[Ln][Ln-1] = -self.Alpha(k, i)
            if(i+1 < self.Nx): # u(x+1,z) inside grid in I
                self.mUt[Ln][Ln+1] = -self.Alpha(k, i)
            if(k-1 >= 0): #u(x,z-1)
                self.mUt[Ln][Ln-self.Nx]= -self.Alpha(k, i)
            if(k+1 < self.Nz): #u(x,z+1)
                self.mUt[Ln][Ln+self.Nx]= -self.Alpha(k, i)

        return self.mUt

    def Independent(self):
        """
        Independent term
        Depends on alpha, velocity and pressure
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
            
            
            self.vId[Ln] = self.Alpha(k, i)*(u0+u1+u2+u3)
            self.vId[Ln] += (2-4*self.Alpha(k, i))*u[1][k][i] - u[0][k][i]

        return self.vId
      
    def SolveNextTime(self):
        """
        Calculate the next time (factorization)
        and update the time stack grids
        """

        try:
            self.mUtfactor
        except :
            self.iter=1
            self.LinearSystem()
            # gets the m factor from the solved system
            self.mUtfactor = ln.lu_factor(self.mUt)

        # modification in the grid due source position
        self.Source(self.iter)
        # As t is in [0, 1, 2] (2nd order)
        # time t in this case is Utime[2]
        # the independent term of the matrix, due the pressure field
        v = self.Independent()

        result = ln.lu_solve(self.mUtfactor, v)
        # reshape the vector to become a matrix again
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
    

    def Source(self, it):
        """
        ( wavelet ) Set the source condition at the pertubation source position.
        ( Si, Sk ) source position
        ( it ) At the given time step. Set t and t-1. (2 order finite diferences)        
        """
        if( it > np.size(self.Wavelet)-1):
            return
           
        self.Utime[1][self.Sk][self.Si] = self.Wavelet[it]        
        ## it should not necessary to replace back the old perturbation value
        ##if(it - 1 < np.size(Wavelet)):
        ##    self.Utime[0][Sk][Si] = Wavelet[it-1]

        return
        
    def Clean():
        """
        De-alloc variables, specially the huge matrix
        used to solve the problem
        """
        self.mUt = None
        self.mUtfactor = None
        gc.collect()
        
        
    def Loop(self):
        """
        Loop through all time steps until (MaxIter)
        saving the matrix at every (Nrec) interactions
        """
        snapiter=0
                
        movie = np.zeros([int(self.MaxIter/self.Nrec), self.Nz, self.Nx])
        # for little problems with the wavelet put initialize as 1
        for i in range(1, self.MaxIter, 1):
            self.SolveNextTime()
            if(i%self.Nrec==0): # every n'th Nrec
                movie[snapiter] = self.Utime[1]
                snapiter+=1
            sys.stdout.write("\r progressing .. %.1f%%" %(100.0*float(i)/self.MaxIter))
            sys.stdout.flush()        
        sys.stdout.write(" done! \n")
                
        return movie
