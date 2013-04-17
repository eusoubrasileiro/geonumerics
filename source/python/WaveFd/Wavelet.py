#!/usr/bin/python
# import filters.py


"""
MAKE SURE THE IPYTHON HAS LOADED THE LAST
VERSION THAT YOU ARE WORKING!!!!
A LOT OF TIME YOU LOST BECAUSE OF THAT!!!!!
"""
#backend = 'gtk'
#import matplotlib
#matplotlib.use(backend)

import sys
sys.path.append('../../python');
# or to not windows guys
# sys.path.append(os.path.join('..','log','python'))
# the above to be able to load the modules bellow

import pylab as py
import numpy as np
from Filters import SincLowPass as _SincWavelet
from Filters import WindowHann

def SincWavelet(N=127, Fc=40, dt=0.0005, plot=False):
    """   
    N is number of samples for the wavelet
    fc is the F maximum frequency
    dt is the sample rate
    Sample the Box Sinc filter simetrically around the zero
    apply hanning window
    """
    print "total wavelet time : %.1f miliseconds" % (dt*N*1000)
    wavelet = _SincWavelet(N, Fc, dt)
    wavelet = wavelet*WindowHann(N)
    #Normalize the maximum amplitude to 1
    wavelet = wavelet/np.max(wavelet)
    return wavelet

def Triangle(N=None, Fc=40.0, Dt=0.001):
    """
    Triangle Wave one Period
    Defined by 1 or 2:
    1) N (half length)
    2) 
        a) Desired Frequency Fc
        b) Sample Rate Dt
    """
    if(N==None):
        N=int(1/float(Fc*Dt))        
    
    t = np.arange(0+1.0/N, 1, 1.0/N)
    y = 1-t
    y = np.append(y, 0.0)
    y_ = 1-t[::-1]
    y_ = np.insert(y_, 0, 0.0)
    
    return np.append(y_, np.append(1, y))


def RickerWavelet(sigma, t):
    wv = 2*np.sqrt(3*sigma)*np.pi**0.25
    return wv* (1-(t/sigma)**2)*np.exp(-0.5*(t/sigma)**2) 


def LinearSin(Fc=40.0, dt=None, plot=False):
    """
    Linear decreasing one period sin(2pi*f) 
    """
    # frequency of niquest limit
    if ( dt > 1/(2.0*Fc) or dt == None):
        dt = 1/(2.0*Fc) 
    
    t = np.arange(0, 1.0/Fc, dt)
    wavelet = np.sin(2.0*np.pi*Fc*t)*(-Fc*t+1.0)
    
    print "total wavelet time : %.1f miliseconds" % (dt*np.size(wavelet)*1000)
    #normalize to 1 the maximum amplitude
    wavelet = wavelet/np.max(wavelet)
    
    # invert the function so, it starts with a small perturbation [::-1]
    return wavelet


def GaussCos(Fc=40.0, dt=None, plot=False):
    """
    TODO
    Linear decreasing one period cos... 
    """    
    t = np.arange(-1.0/Fc, 0, dt)
    wavelet = np.cos(2.0*np.pi*Fc*t)*(Fc*t+1.0)
    t = np.arange(0, 1.0/Fc, dt)
    wavelet = np.append(wavelet, np.cos(2.0*np.pi*Fc*t)*(-Fc*t+1.0) )


class SourceWavelet:
    """
    Represents the source wavelet, characteristic of the source
    creating the perturbation wave field
    """
    def __init__(self, Fc=40.0, Type=None):
        """
        Fc - central frequency
        
        Type can be:
        LinearSin - 1) Linear decreasing one period sin(2pi*f) 
        """
        
        self.Fc = Fc

        if ( Type == None):            
                self._Type = SincWavelet
        else:
            if ( Type == None):
                self._Type = LinearSin
        
    def Samples(self, dt):
        """
        get the values itself based on the predefined Fc
        and the passed dt
        """

        if( self._Type == SincWavelet):
            s = self._Type(1/(self.Fc*dt), Fc=self.Fc, dt=dt)
        else:
            s = self._Type(Fc=self.Fc,dt=dt)

        return s

