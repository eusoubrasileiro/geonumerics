import numpy as np
from BaseWaveField import BaseWave2DField
import pylab as py

"""
Tenth tentative using polynomial interpolations
to approximate the derivatives
now uses the 0 padding before calculating the derivatives.
"""

def FourierDerivative(f):
    """
    this derivatie just works for periodic 2*pi multiple series
    have to figure out how to make that work for any function
    """
    N = np.size(f)
    n = np.arange(0,N)
    # df discrete differential operator
    df = np.complex(0,1)*py.fftshift(n-N/2)
    dfdt = py.ifft( df*py.fft(f) )  
    return py.real(dfdt)
    


# amazing aproximation but very timing consuming
# and its also limited to the number of points that you use
def ChebyshevDerivative(x, y, m=4):
    """
    Fits a Chebyshev poly of degree (deg)
    in (y) and calculate the 1st derivative.
    
    (2*m+1) points are used to calculate the derivative
    in the borders previous/after points are used.
    """
    deg = 2*m # we always need 2*m+1 points for a degree 2*m
    n = np.size(x)
    if(n != np.size(y)
        or n < deg 
        or n < m*2+1):
        raise Exception(" size/degree/window inconsistence  ")

    # Chebyshev polys in the interval
    npoly = n-m*2
    Cheb_Coef = np.zeros([npoly, deg+1])          
    for i in range(npoly):
        Cheb_Coef[i] = cheby.chebfit(x[i:i+2*m+1], y[i:i+2*m+1], deg)

    # der coeficients    
    Cheb_Coefder = np.zeros([npoly, deg])          
    for i in range(npoly):
        Cheb_Coefder[i] = cheby.chebder(Cheb_Coef[i], 1)
   
    dydx = np.zeros(n)    
    dydx[0:m] = cheby.chebval(x[0:m], Cheb_Coefder[0]) #begin
    dydx[-m:] = cheby.chebval(x[-m:], Cheb_Coefder[npoly-1]) #end    
    for i in range(m, n-m): #middle
        dydx[i] = cheby.chebval(x[i], Cheb_Coefder[i-m])  

    return dydx

def SplineDerivative(x, y, k=3):
    """
    Fits a spline of degree (k)
    in (y) and calculate the 1st derivative.
    """
    n = np.size(x)
    if(n != np.size(y)
        or n < k):
        raise Exception(" size/degree inconsistence  ")

    Spline_Interp = interp.InterpolatedUnivariateSpline(x, y, k=k)
    dydx = np.zeros(np.size(x))
    # first derivative [1]
    for i in range(np.size(x)):
        dydx[i] = Spline_Interp.derivatives(x[i])[1]

    return dydx

def SplineDerivativeZero(Ds, y, k=3):
    """
    padd k zeros in the borders and calculate the derivative
    """
    # padd with zeros before and after size k
    padd = k+3
    y_ = np.append(np.zeros(padd), np.append(y, np.zeros(padd)))
    x = np.arange(0.0, Ds*np.size(y_), Ds)
    Spline_Interp = interp.InterpolatedUnivariateSpline(x, y_, k=k)
    dydx = np.zeros(np.size(y_))
    # first derivative [1]
    for i in range(np.size(x)):
        dydx[i] = Spline_Interp.derivatives(x[i])[1]

    return (dydx[padd:])[:-padd]

class Exp2DFourierWave(BaseWave2DField):
    r"""
    Explicit 2D wave equation with 2 order centered in time
    and using spline for derivatives in space        
    """

    def __init__(self,
                nx, 
                nz,
                ds,
                dt,
                velocity,
                rho,
                sx,
                sz,
                maxiter,
                nrec=1,
                wavelet=None):
        r"""
        Initialize a new wave equation field explicit centered differences time
        and spline for space.
        
        * nx       : number of discretization in x
        * nz       : number of discretization in z
        * ds       : dx=dz=ds grid spacing 
        * dt       : time step - e.g. seconds
        * velocity : 2d velocity distribution
        * rho      : 2d density distribution
        * sx/sz    : source wavelet position in indexes (i, k)
        * maxiter  : total iterations
        * nrec     : recording interval 1 equals time step 
        * wavelet  : source wavelet function applied at the position (sx, sz)
          must have sample rate equal to dt
        """
        # creates the constant density part first
        super(Exp2DFourierWave, self).__init__(nx, nz, ds, dt, velocity, sx, sz, maxiter, nrec, wavelet)
        
        # centered differences in time, add previous time
        self.Uprevious = np.zeros([self.Nz, self.Nx])
        # Rho field for each (x, z) point
        # Rho doesnt vary with time
        self.Rho = np.zeros([self.Nz, self.Nx])
        # if a matrix of velocity is passed fills it        
        # setts up the velocity field
        if(type(rho) is np.ndarray):
            if(np.shape(rho) == (self.Nz, self.Nx) ):
                self.Rho = rho
        # if not put a constant velocity
        else:
            self.Rho[:][:] = rho


    def SolveNextTime(self):

        try:
            self.tstep += 1
        except :
            self.tstep = 0

        self.Source(self.tstep)
        
        # calculate derivatives using fourier
        # derivatives in x and y direction
        d2dx = np.zeros([self.Nz, self.Nx]) 
        d2dz = np.zeros([self.Nz, self.Nx])

        # xline by xline
        for i in range(self.Nz):        
            d2dx[i] = FourierDerivative(self.Ucurrent[i]) # first derivative
            d2dx[i] = FourierDerivative(d2dx[i]/self.Rho[i]) # second with Rho

        # to easy calculate on the z direction, transpose the matrixes
        self.Ucurrent = self.Ucurrent.transpose()
        self.Rho = self.Rho.transpose()
        d2dz = d2dz.transpose()
        # zline by zline
        for i in range(self.Nx):        
            d2dz[i] = FourierDerivative(self.Ucurrent[i]) # first derivative
            d2dz[i] = FourierDerivative(d2dz[i]/self.Rho[i]) # second with Rho
        # transpose back
        self.Ucurrent = self.Ucurrent.transpose()
        self.Rho = self.Rho.transpose()
        d2dz = d2dz.transpose()

        # derivatives combination + source
        LP = d2dx + d2dz      
        # simple centered differences in time 2nd order
        # matrix operations
        self.Ufuture = ((self.Vel*self.Dt)**2)*LP*self.Rho+2*self.Ucurrent-self.Uprevious     

        # update stack times
        self.Uprevious = self.Ucurrent
        self.Ucurrent = self.Ufuture        
        
        return self

    
    
