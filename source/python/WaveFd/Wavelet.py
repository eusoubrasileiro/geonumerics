r"""

**Collection of source functions or so called Wavelets**

.. note:: 
    Remember that nyquist frequency is :math:`\frac{1}{2dt}` where
    `dt` is sample rate.

.. note::
    All functions *return* the sampled wavelet

.. todo::

    Done, wavelets should be sample at run time, not prior. Using __call__ solves
    most of the problem while also the wavelet the can directly used in the analytical
    solution, used to avoid big gradients. 

"""
import sys
sys.path.append('../../python');
import numpy as np

class RickerSource(object):
    r"""
    Ricker Wavelet:
    :math:`A = (1-2 \pi^2 f^2 t^2) e^{-\pi^2 f^2 t^2}`

    * fc        : maximum desired frequency
    * amp       : source strength 
    * delay     : to turn the source causal recommend 1/fc (starting at zero)
    """

    def __init__(self, fc, amp, delay=0):
        self.fc = fc
        self.amp = amp
        self.delay = delay

    def __call__(self, t):
        t = t-self.delay
        cx = (np.pi*self.fc*t)**2
        return  self.amp*(1-2*cx)*np.exp(-cx) 



from Filters import SincLowPass, WindowHann, FilterSize

# default relative bandwidth transition for defining wavelet size
RBTW = 0.40


def Sinc(fc, dt, n=None):
    r"""   
    Sinc Wavelet source.
    Sample a Sinc Box filter simetrically around the zero
    applying a hanning window

    * fc : maximum frequency
    * dt : sample rate
    * n  : number of samples for the wavelet (odd)
    """
    if(n==None):
        n=FilterSize(fc, dt, RBTW)
    wavelet = SincLowPass(n, fc, dt)
    wavelet = wavelet*WindowHann(n)
    #Normalize the maximum amplitude to 1
    wavelet = wavelet/np.max(wavelet)
    return wavelet

def Triangle(fc, dt, n=None):
    r"""
    Triangle Wave one Period.
    Defined by frequency and sample rate or by size

    * fc        : maximum desired frequency
    * dt        : sample rate
    * n         : half length of triangle    
    """
    if(n==None):
        n=int(1/float(fc*dt))

    t = np.arange(0+1.0/n, 1, 1.0/n)
    y = 1-t
    y = np.append(y, 0.0)
    y_ = 1-t[::-1]
    y_ = np.insert(y_, 0, 0.0)
    
    return np.append(y_, np.append(1, y))



def Ricker(fc, dt, n=None):
    r"""
    Ricker Wavelet:
    :math:`A = (1-2 \pi^2 f^2 t^2) e^{-\pi^2 f^2 t^2}`

    * fc        : maximum desired frequency
    * dt        : sample rate
    """
    # n is specified based on closer to zero 
    # 0.0 = (1-2 pi^2 fc^2 t^2) e^(-pi^2 fc^2 t^2)
    # that gives {t = -/+ 1/(sqrt(2) pi fc)}
    
    if(n==None):
        n=FilterSize(fc, dt, RBTW)
    t = np.arange(-dt*(n-1)/2,(dt*(n-1)/2)+dt, dt)
    # avoid division by zero
    #t[np.size(t)/2]=np.NaN
    ricker = (1-2*(np.pi*fc*t)**2)*np.exp(-(np.pi*fc*t)**2) 
    return  ricker/np.max(ricker)


def LinearSin(fc, dt, n=None):
    r"""
    Linear decreasing one period sin(2pi*f) 
    
    * fc : maximum desired frequency
    * dt : sample rate
    """
    # frequency of niquest limit
    if ( dt > 1/(2.0*fc)):
        dt = 1/(2.0*fc) 
    
    t = np.arange(0, 1.0/fc, dt)
    wavelet = np.sin(2.0*np.pi*fc*t)*(-fc*t+1.0)
    #normalize to 1 the maximum amplitude
    wavelet = wavelet/np.max(wavelet)
    
    # invert the function so, it starts with a small perturbation [::-1]
    return wavelet





