import numpy as np
from BaseWaveField import BaseWave2DField
from WaveFd.Utils import centered2ndfiniteweight4th
# pyx optimization
import pyximport; pyximport.install()
import _cExp2DWave

class Exp2DWave(BaseWave2DField):
    r"""
    Simple explicit wave equation constant density accoustic, 
    Finite differences:

    * 4th error order centered in space
    * 1st error order backward in time

    .. todo::
    
    * Very Slow : try using cython for loop said to be 8x faster!!!
    Something is wrong in the implementation, the convergence criteria is
    correc by Jing-Bo Chen or something else is wrong.
    
    Convergence
     
    .. math::
     
         \Delta t \leq \frac{2 \Delta s}{ V \sqrt{\sum_{a=-N}^{N} (|w_a^1| + |w_a^2|)}}
    
    Where w_a are the centered differences weights
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
                wavelet=None):
        r"""
        Initialize a new wave equation field explicit centered differences time

        * nx       : number of discretization in x
        * nz       : number of discretization in z
        * ds       : dx=dz=ds grid spacing 
        * dt       : time step - e.g. seconds
        * velocity : 2d velocity distribution
        * sx/sz    : source wavelet position in indexes (i, k)
        * maxiter  : total iterations
        * nrec     : recording interval 1 equals time step 
        * wavelet  : source wavelet function applied at the position (sx, sz)
          must have sample rate equal to dt
        """
        # creates the constant density part first
        super(Exp2DWave, self).__init__(nx, nz, ds, dt, velocity, sx, sz, maxiter, nrec, wavelet)
        # backward/forward differences in time, add previous time
        self.Uprevious = np.zeros([self.Nz, self.Nx])
        
        neededt = self.Stability()
        
        if(dt > neededt):
            print "for stability the time step should be smaller than", neededt
            return 

    def Stability(self):
        r"""
        Using Von neuman stability analysis
        
        .. math:
     
         \Delta t \leq \frac{2 \Delta s}{ V \sqrt{\sum_{a=-N}^{N} (|w_a^1| + |w_a^2|)}}
    
        Where w_a are the centered differences weights. And V is the maximum V. 
        
        * Returns maximum value allowed for \Delta t
        """
        
        vmax = 0.0
        
        for k in range(self.Nz):
            for j in range(self.Nx):
                if ( vmax < self.Vel[k][j]):
                    vmax =  self.Vel[k][j]
        
        sumweights = np.abs(centered2ndfiniteweight4th).sum()
        
        return 2*self.Ds/(vmax*np.sqrt(sumweights))

    # def NextTime(self, k, i):

    #     u = self.Ucurrent
    #     # u0k u1k*uik*u3k u4k     
    #     # Boundary fixed 0 outside        
    #     u0k=u1k=u3k=u4k=0.0
    #     ui0=ui1=ui3=ui4=0.0
    #     uik = u[k][i]      
          
    #     if(i-2 > -1):
    #         u0k = u[k][i-2]
    #     if(i-1 > -1):
    #         u1k = u[k][i-1]
    #     if(i+1 < self.Nx):
    #         u3k = u[k][i+1]            
    #     if(i+2 < self.Nx):
    #         u4k = u[k][i+2]
    #     if(k-2 > -1):
    #         ui0 = u[k-2][i]
    #     if(k-1 > -1):
    #         ui1 = u[k-1][i]
    #     if(k+1 < self.Nz):
    #         ui3 = u[k+1][i]            
    #     if(k+2 < self.Nz):
    #         ui4 = u[k+2][i]

    #     d2u_dx2 = (-u0k+16*u1k-30*uik+16*u3k-u4k)/12.0
    #     d2u_dz2 = (-ui0+16*ui1-30*uik+16*ui3-ui4)/12.0
    #     Uikfuture = (d2u_dx2+d2u_dz2)*(self.Dt*self.Vel[k][i]/self.Ds)**2
    #     Uikfuture += 2*self.Ucurrent[k][i]-self.Uprevious[k][i]
         
    #     return Uikfuture


    # def SolveNextTime(self):

    #     try:
    #             self.tstep += 1
    #     except :
    #             self.tstep = 0
        
    #     self.Source(self.tstep)            

    #     for k in range(self.Nz):
    #         for j in range(self.Nx):
    #             self.Ufuture[k][j] = self.NextTime(k, j)

    #     # make the update in the time stack
    #     self.Uprevious = self.Ucurrent
    #     self.Ucurrent = self.Ufuture
    
    #     return self.Ufuture


    def SolveNextTime(self):

        try:
                self.tstep += 1
        except :
                self.tstep = 0
        
        self.Source(self.tstep)   

        _cExp2DWave.SolveUfuture(self.Ufuture, self.Ucurrent, self.Uprevious, self.Vel, self.Nz, self.Nx, self.Ds, self.Dt)

        # make the update in the time stack
        self.Uprevious = self.Ucurrent
        self.Ucurrent = self.Ufuture
    
        return self.Ufuture

#  CYTHON!
