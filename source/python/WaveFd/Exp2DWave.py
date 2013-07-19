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
    correct by Jing-Bo Chen. Problem is big gradients causing instability due
    source?? 
    
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
                wavelet,
                nrec=1):
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
        """

        # creates the constant density part first
        super(Exp2DWave, self).__init__(nx, nz, ds, dt, velocity, sx, sz, maxiter, wavelet, nrec)
        # backward/forward differences in time, add previous time
        self.Uprevious = np.zeros([self.Nz, self.Nx])
        
        neededt = self._MaxDt()
        
        if(dt > neededt):
            print "for stability the time step should be smaller than", neededt
            return 
        
        self._Summary()
        
    def _NiquestSpatial(self):
        r"""
        Spatial niquest in wave number based on spatial sample rate
        """
        return 1.0/(2*self.Ds)

    def _SpatialAlias(self):
        r"""
        Returns True or False
        False if the maximum wavelength is bigger than Niquest so acceptable
        True otherwise
        """
        print ((self.vmin/self.Wavelet.fc) < 1.0/self._NiquestSpatial())

    def _GridPointsbyWavelenght(self):
        r"""
        Alford et all
        measures the number of points (based on spatial sample rate)
        by wave length. Wavelenght calculated from source and velocity
        wlenght = velocity/frequency
        """
        return (self.vmin/self.Wavelet.fc)/self.Ds 

    def _Summary(self):
        r"""
        Summary all information for this simulation
        """
        print "R :", self.vmax*self.Dt/self.Ds, " of allowed < : ", np.sqrt(3)/np.sqrt(8)
        print "Dt :", self.Dt, " of allowed < : ", self._MaxDt()
        print "Points by wavelenght: ", self._GridPointsbyWavelenght(), " recommended > 5"
        print "Is there spatial Alias (based on frequency)? ", self._SpatialAlias()
        
        return

    def _MaxDt(self):
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
                if ( vmax < self.Vel[k, j]):
                    vmax =  self.Vel[k, j]
        
        sumweights = np.abs(centered2ndfiniteweight4th).sum()
        
        return 2*self.Ds/(vmax*np.sqrt(sumweights))

    def SolveNextTime(self):

        try:
                self.tstep += 1
        except :
                self.tstep = 0
        
        _cExp2DWave.SolveUfuture(self.Ufuture, self.Ucurrent, self.Uprevious, self.Vel, self.Nz, self.Nx, self.Ds, self.Dt)

        # source position center grid
        # smooth region around the center of the grid +3-3
        # exact analytical solution circular
        # 2.5 wavelength is enough for the source
        if(self.tstep*self.Dt < 2.5/self.Wavelet.fc):
            dsm = 2
            c = self.Vel[self.Sk][self.Si] # uses central velocity for exact solution
            R2 = (c*self.Dt)**2
            for nz in range(max(self.Sk-dsm,0),min(self.Sk+dsm+1,self.Nz)):
                for nx in range(max(self.Si-dsm,0),min(self.Si+dsm+1,self.Nx)):
                    dz = nz - self.Sk 
                    dx = nx - self.Si
                    t = self.tstep*self.Dt
                    r = np.sqrt(dz**2+dx**2)
                    if(nx == self.Si and nz == self.Sk):
                        self.Ufuture[nz,nx] -= R2*self.Wavelet(t)
                    else:        
                        self.Ufuture[nz,nx] -= R2*self.Wavelet(t-r/c)/(2*np.pi*r)
            
        # make the update in the time stack
        self.Uprevious[:][:] = self.Ucurrent[:][:]
        self.Ucurrent[:][:] = self.Ufuture[:][:]
    
        return self.Ufuture

#  CYTHON!
