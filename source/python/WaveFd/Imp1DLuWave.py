import numpy as np
import scipy.linalg as ln
from BaseWaveField import BaseWave1DField


class Imp1DLuWave(BaseWave1DField):
    r"""
    Simple implicit 1D wave equation LU decomposition wave equation constant
    density, Finite differences:

    * 2nd order centered in space
    * 1st order backward in time
    """
    def __init__(self,
                nx, 
                ds,
                dt,
                velocity,
                sx,
                maxiter,
                nrec=1,
                wavelet=None):
        r"""
        Initialize a new wave equation field

        * nx       : number of discretization in x
        * ds       : dx=ds grid spacing 
        * dt       : time step - e.g. seconds
        * velocity : 1d velocity distribution
        * sx       : source wavelet position in indexes i
        * maxiter  : total iterations
        * nrec     : recording interval 1 equals time step 
        * wavelet  : source wavelet function applied at the position sx
          must have sample rate equal to dt
        """
        super(Imp1DLuWave, self).__init__(nx, ds, dt, velocity, sx, maxiter, nrec, wavelet)
        # backward/forward differences in time, add previous time
        self.Uprevious = np.zeros(self.Nx)


    def LinearSystem(self):
        r"""
        Assembly linear system
        """
        
        # assembly matrix of linear system
        # to solve u(t) based on u(t-1) and u(t-2)
        # the matrix includes all future values of u
        # in the entire grid, so size is the number of cells
        # start with zeros that is also the countour condition u(t)=0
        self.mUt = np.zeros([self.Nx, self.Nx])
        # assembly linear system, the linear system
        # ignores external part of the grid = locked boundary
        # ln go through all the cells in the grid Ut
        # each cell gives one equation (line)
        for Ln in range(0, self.Nx): 
            # 1.0*u(x-1) + gama(x)*u(x) + 1.0*u(x+1) 
            # turn the indices to the one of original matrix
            i = Ln
            # just major diagonal
            self.mUt[i][i] = -2.0-(self.Ds/(self.Dt*self.Vel[i]))**2

            if(i-1 >= 0): # u(x-1) inside grid in I
                self.mUt[Ln][Ln-1] = 1.0
            if(i+1 < self.Nx): # u(x+1) inside grid in I
                self.mUt[Ln][Ln+1] = 1.0
            
        return self.mUt


    def Independent(self):
        r"""
        Independent term
        """
        #independent term, where the previous times goes in
        self.vId = np.zeros([self.Nx])
        # fill the independent vector
        for Ln in range(0, self.Nx): 
            # turn the indices to the one of original matrix            
            # -2 u(x,t-1) + u(x,t-2)
            self.vId[Ln] = -2*self.Ucurrent[Ln]+self.Uprevious[Ln]
            # / (Ds/(Dt*Vel(x,z)))**2
            self.vId[Ln] *= (self.Ds/(self.Dt*self.Vel[Ln]))**2

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

        self.Source(self.tstep)
        # in time
        # As t is in [0, 1, 2] (2nd order)
        # time t in this case is Utime[2]
        v = self.Independent()
        result = ln.lu_solve(self.mUtfactor, v)
        # reshape the vector to became a matrix again
        self.Ufuture = result
        # make the update in the time stack
        # before [t-2, t-1,  t]
        # after  [t-1, t,  t+1]
        # so t-2 receive t-1 and etc.
        self.Uprevious = self.Ucurrent 
        self.Ucurrent = self.Ufuture
        
        return self.Ufuture