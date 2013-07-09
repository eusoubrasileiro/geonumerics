"""
Implicit 2D wave equation 

"""

import numpy as np
import scipy.linalg as ln
from BaseWaveField import BaseWave2DField

class Imp2DLuWave(BaseWave2DField):
    r"""
    Implicit wave equation constant density accoustic, 
    Finite differences:
    
    * 2nd error order centered in space
    * 1st error order backward in time
    
    .. note:

        Using Lu decomposition and factor

    .. note:
    
        Cholesky_band requires LinAlgError: 1-th leading minor not positive definite
        cho_factor gives wrong results as time varies
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
                 wavelet,
                 nrec=1
                 ):
        r"""
        Initialize a new wave equation field implicit centered differences second order

        * nx       : number of discretization in x
        * nz       : number of discretization in z
        * ds       : dx=dz=ds grid spacing
        * dt       : time step - e.g. seconds
        * velocity : 2d velocity distribution
        * sx/sz    : source wavelet position in indexes (i, k)
        * maxiter  : total iterations
        * wavelet  : source wavelet function applied at the position (sx, sz)
        * nrec     : recording interval 1 equals time step 
        """
        super(Imp2DLuWave, self).__init__(nx, nz, ds, dt, velocity, sx, sz, maxiter, wavelet, nrec)
        self.Uprevious = np.zeros([self.Nz, self.Nx])


    def Gamma(self, k, i):
        r"""
        .. math:

        \gamma = -4 -r^{-2}
        """
        return -(4 +self.R_(k, i)**2)

    def R(self, k, i):
        r"""
        .. math:
        
         r = \frac{\Delta t  V_{jk}^n}{ \Delta s}        
        """
        return self.Dt * self.Vel[k][i]/ self.Ds

    def R_(self, k, i):
        r"""
        .. math:
        
         r = \frac{\Delta s}{\Delta t  V_{jk}^n}        
        """
        return self.Ds /(self.Dt * self.Vel[k][i]) 

    def LinearSystem(self):
        r"""
        Assembly linear system
        Depends on Velocity field and Gamma
        """
        # assembly matrix of linear system
        # to solve u(t) based on u(t-1) and u(t-2)
        # the matrix includes all future values of u
        # in the entire grid, so size is the number of cells
        # start with zeros that is also the boundary condition u(t)=0
        self.mUt = np.zeros([self.Nz*self.Nx, self.Nz*self.Nx])

        # assembly linear system, the linear system
        # ignores external part of the grid = locked boundary
        # ln go through all the cells in the grid Ut
        # each cell gives one equation (line)
        for Ln in range(0, self.Nz*self.Nx, 1):
            # 1.0*u(x-1,z) + Gamma(x,z)*u(x,z) + 1.0*u(x+1,z) + 1.0*u(x,z-1) + 1.0*u(x,z+1)
            # turn the indices to the one of original matrix
            i = Ln%self.Nx
            k = Ln/self.Nx

            self.mUt[Ln][Ln] = self.Gamma(k, i)
            #is this right?
            if(i-1 >= 0): # u(x-1,z) inside grid in I
                self.mUt[Ln][Ln-1] = 1.0
            if(i+1 < self.Nx): # u(x+1,z) inside grid in I
                self.mUt[Ln][Ln+1] = 1.0
            if(k-1 >= 0): #u(x,z-1)
                self.mUt[Ln][Ln-self.Nx]= 1.0
            if(k+1 < self.Nz): #u(x,z+1)
                self.mUt[Ln][Ln+self.Nx]= 1.0

        return self.mUt


    def Independent(self):
        r"""
        Independent term
        Depends on Gamma, velocity and pressure
        """
        #independent term, where the previous times goes in
        self.vId = np.zeros([self.Nz*self.Nx])

        # fill the independent vector
        for Ln in range(0, self.Nz*self.Nx, 1):
            # turn the indices to the one of original matrix
            i = Ln%self.Nx
            k = Ln/self.Nx
            # boundary locked
            self.vId[Ln] = (self.Uprevious[k][i]-2*self.Ucurrent[k][i])*(self.R_(k, i)**2)

        return self.vId

    def SolveNextTime(self):
        r"""
        Calculate the next time (factorization)
        and update the time stack grids
        """

        try:
            self.tstep += 1
        except :
            self.tstep = 0
            self.LinearSystem()
            # gets the m factor from the solved system
            self.mUtfactor = ln.lu_factor(self.mUt)

        # modification in the grid due source position
        self.Ucurrent[self.Sk][self.Si] = self.Wavelet(tstep*self.Dt)   
        # As t is in [0, 1, 2] (2nd order)
        # time t in this case is Utime[2]
        # the independent term of the matrix, due the pressure field
        v = self.Independent()

        result = ln.lu_solve(self.mUtfactor, v)
        # reshape the vector to become a matrix again
        self.Ufuture = np.reshape(result, (self.Nz, self.Nx))

        # make the update in the time stack
        # before [t-2, t-1,  t]
        # after  [t-1, t,  t+1]
        # so t-2 receive t-1 and etc.
        # make the update in the time stack
        self.Uprevious = self.Ucurrent
        self.Ucurrent = self.Ufuture        
        
        return self.Ufuture


from pysparse.sparse import spmatrix
from pysparse.direct import umfpack


class Imp2DLuSparseWave(Imp2DLuWave):
    r"""
    Implicit wave equation constant density accoustic, 
    Finite differences:
    
    * 2nd error order centered in space
    * 1st error order backward in time
    
    .. note:

        Using Pysparse Lu decomposition and factor. Huge memory drop.
        Eg. for a 200x300 grid, memory usage for the linear system is 0.0083%
        compared with the Imp2DLuWave matrix above. Or ~12*10^3 smaller  

    """

    def LinearSystem(self):
        r"""
        Assembly linear system
        Depends on Velocity field and Gamma
        """
        # assembly matrix of linear system
        # using pysparse optimized matrix non zero elements 5*M         
        self.mUt = spmatrix.ll_mat(self.Nz*self.Nx, self.Nz*self.Nx, 5*self.Nz*self.Nx-2*self.Nz-2*self.Nx)

        for Ln in range(0, self.Nz*self.Nx, 1):
            # 1.0*u(x-1,z) + Gamma(x,z)*u(x,z) + 1.0*u(x+1,z) + 1.0*u(x,z-1) + 1.0*u(x,z+1)
            # turn the indices to the one of original matrix
            i = Ln%self.Nx
            k = Ln/self.Nx

            self.mUt[Ln,Ln] = self.Gamma(k, i)
            #is this right?
            if(i-1 >= 0): # u(x-1,z) inside grid in I
                self.mUt[Ln,Ln-1] = 1.0
            if(i+1 < self.Nx): # u(x+1,z) inside grid in I
                self.mUt[Ln,Ln+1] = 1.0
            if(k-1 >= 0): #u(x,z-1)
                self.mUt[Ln,Ln-self.Nx]= 1.0
            if(k+1 < self.Nz): #u(x,z+1)
                self.mUt[Ln,Ln+self.Nx]= 1.0

        return self.mUt


    def SolveNextTime(self):
        r"""
        Calculate the next time (factorization)
        and update the time stack grids
        """

        try:
            self.tstep += 1
        except :
            self.tstep = 0
            self.LinearSystem()
            self.mUtLU = umfpack.factorize(self.mUt, strategy="UMFPACK_STRATEGY_SYMMETRIC")
            # gets the m factor from the solved system
        # modification in the grid due source position
        self.Ucurrent[self.Sk][self.Si] = self.Wavelet(tstep*self.Dt)   
        # As t is in [0, 1, 2] (2nd order)
        # time t in this case is Utime[2]
        # the independent term of the matrix, due the pressure field
        v = self.Independent()
        result = np.empty(self.Nx*self.Nz)
        self.mUtLU.solve(v, result)
        # reshape the vector to become a matrix again
        self.Ufuture = np.reshape(result, (self.Nz, self.Nx))

        # make the update in the time stack
        # before [t-2, t-1,  t]
        # after  [t-1, t,  t+1]
        # so t-2 receive t-1 and etc.
        # make the update in the time stack
        self.Uprevious = self.Ucurrent
        self.Ucurrent = self.Ufuture        
        
        return self.Ufuture