#!/usr/bin/python

"""

Collection of source functions or so called Wavelets.
Remember that nyquist frequency is 1/(2 sample rate)

"""
import sys
sys.path.append('../../python');

import numpy as np
from Filters import SincLowPass
from Filters import WindowHann

def SincWavelet(n=127, fc=40, dt=0.0005):
    """   
    Sinc Wavelet source.
    Sample a Sinc Box filter simetrically around the zero
    applying a hanning window
    
    n  : number of samples for the wavelet
    fc : maximum frequency
    dt : sample rate
    """
    print "total wavelet time : %.1f miliseconds" % (dt*n*1000)
    wavelet = SincLowPass(n, fc, dt)
    wavelet = wavelet*WindowHann(n)
    #Normalize the maximum amplitude to 1
    wavelet = wavelet/np.max(wavelet)
    return wavelet

def Triangle(n=None, fc=40.0, dt=0.001):
    """
    Triangle Wave one Period.
    Defined by size or by frequency and sample rate
    
    n  : half length of triangle    
    fc : maximum desired frequency
    dt : sample rate
    """
    if(n==None):
        n=int(1/float(fc*dt))        
    
    t = np.arange(0+1.0/n, 1, 1.0/n)
    y = 1-t
    y = np.append(y, 0.0)
    y_ = 1-t[::-1]
    y_ = np.insert(y_, 0, 0.0)
    
    return np.append(y_, np.append(1, y))


#TODO: finish
# def RickerWavelet(fc, dt):
#     """
#     Ricker Wavelet    
#     A = (1-2 \pi^2 f^2 t^2) e^{-\pi^2 f^2 t^2}
#     Side lobes :
#     \pm \frac{\sqrt{3/2}}{f\pi}
#      
#     fc : maximum desired frequency
#     dt : sample rate
#     """
#     
#     n=int(1/float(fc*dt))
#     period=1.0/fc    
#     t = np.arange(-n, 1, 1.0/n)
#     ricker = (1-2*(np.pi*fc*t)**2)*np.exp(-(np.pi*fc*t)**2) 
#     return  ricker/np.max(ricker)


def LinearSin(fc=40.0, dt=None):
    """
    Linear decreasing one period sin(2pi*f) 
    
    fc : maximum desired frequency
    dt : sample rate
    """
    # frequency of niquest limit
    if ( dt > 1/(2.0*fc) or dt == None):
        dt = 1/(2.0*fc) 
    
    t = np.arange(0, 1.0/fc, dt)
    wavelet = np.sin(2.0*np.pi*fc*t)*(-fc*t+1.0)
    
    print "total wavelet time : %.1f miliseconds" % (dt*np.size(wavelet)*1000)
    #normalize to 1 the maximum amplitude
    wavelet = wavelet/np.max(wavelet)
    
    # invert the function so, it starts with a small perturbation [::-1]
    return wavelet


# def GaussCos(Fc=40.0, dt=None, plot=False):
#     """
#     TODO
#     Linear decreasing one period cos... 
#     """    
#     t = np.arange(-1.0/Fc, 0, dt)
#     wavelet = np.cos(2.0*np.pi*Fc*t)*(Fc*t+1.0)
#     t = np.arange(0, 1.0/Fc, dt)
#     wavelet = np.append(wavelet, np.cos(2.0*np.pi*Fc*t)*(-Fc*t+1.0) )
