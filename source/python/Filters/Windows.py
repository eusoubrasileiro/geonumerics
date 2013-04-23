'''
Created on Apr 23, 2013

@author: andre
'''
import numpy

"""

Note for all windows size n must be odd

"""

def WindowHann(n):
    """
    Hanning window working for smooth the frequency response
    """
    return numpy.sin(numpy.arange(0, n, 1)*numpy.pi/(n-1));

# seams that those two bellow are broken
# for that try N=5, 7 and the results will be weird

def WindowGaussHann(n, sg=0.25):
    """
    GaussianHanning Window working for smooth the frequency response
    sg < 0.5
    """
    # since Gaussian extends to the infinity it must be cutted
    hw = WindowHann(n);
    return hw*numpy.exp(-0.5*((numpy.arange(0, n, 1)-(n-1)/2)/(sg*(n-1)/2))**2)
    
    
def WindowTukey(n, aph=0.20):
    """
    when tukey aph value = zero it becames a rectangular window
    when tukey aph value = 1 it becames a Hann window
    size of the plateu is (N-1)(1-aph)
    """ 
    nbg = int(aph*(n-1.)/2.)
    nmd = int((n-1.)*(1.-aph))

    bg = 0.5*(1.+numpy.cos(numpy.pi*(-1.+(2.*numpy.arange(nbg)/(aph*(n-1.))))))
    md = 1 + numpy.zeros(nmd)
    en = 0.5*(1.+numpy.cos(numpy.pi*(1.-2./aph+(2.*numpy.arange(nbg+nmd,n)/(aph*(n-1))))))

    # for simplicity and because I dont want to work now...
    return numpy.append(bg, numpy.append(md, en))