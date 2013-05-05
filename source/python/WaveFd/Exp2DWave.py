import numpy as np
from BaseWaveField import BaseWave2DField

class Exp2DWave(BaseWave2DField):
    r"""
    Simple explicit wave equation constant density accoustic, 
    Finite differences:

    * 4th error order centered in space
    * 1st error order backward in time

    .. todo::

     * Convergence??Dt must be very small?? how small?? to calculate!!
     * Very Slow : try using cython for loop said to be 8x faster!!!

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
                nrec=5,
                wavelet=None):
        r"""
        Initialize a new wave equation field explicit centered differences time
        and spline for space.

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
        

    def NextTime(self, i, k):

        u = self.Ucurrent
        # u0k u1k*uik*u3k u4k     
        # Boundary fixed 0 outside        
        u0k=u1k=u3k=u4k=0.0
        ui0=ui1=ui3=ui4=0.0
        uik = u[k][i]      
          
        if(i-2 > -1):
            u0k = u[k][i-2]
        if(i-1 > -1):
            u1k = u[k][i-1]
        if(i+1 < self.Nx):
            u3k = u[k][i+1]            
        if(i+2 < self.Nx):
            u4k = u[k][i+2]
        if(k-2 > -1):
            ui0 = u[k-2][i]
        if(k-1 > -1):
            ui1 = u[k-1][i]
        if(k+1 < self.Nz):
            ui3 = u[k+1][i]            
        if(k+2 < self.Nz):
            ui4 = u[k+2][i]

        d2u_dx2 = (-u0k+16*u1k-30*uik+16*u3k-u4k)/12.0
        d2u_dz2 = (-ui0+16*ui1-30*uik+16*ui3-ui4)/12.0
        Uikfuture = (d2u_dx2+d2u_dz2)*(self.Dt*self.Vel[k][i]/self.Ds)**2
        Uikfuture += 2*self.Ucurrent[k][i]-self.Uprevious[k][i]
         
        return Uikfuture


    def SolveNextTime(self):

        try:
                self.tstep += 1
        except :
                self.tstep = 0
        
        self.Source(self.tstep)            

        for k in range(self.Nz):
            for i in range(self.Nx):
                self.Ufuture[k][i]=self.NextTime(k,i)

        # make the update in the time stack
        self.Uprevious = self.Ucurrent
        self.Ucurrent = self.Ufuture
    
        return self.Ufuture
