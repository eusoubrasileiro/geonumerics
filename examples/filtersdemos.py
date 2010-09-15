#-------------------------------------------------------------------------------
# Name:        filters demos
# Purpose:     some examples using the filters.py library
#
# Author:      andre
#
# Created:     21/08/2010
# Copyright:   (c) andre 2010
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

import sys
sys.path.append('..\\source\\python');
# or to not windows guys
# sys.path.append(os.path.join('..','log','python'))
# the above to be able to load the modules bellow
import samplesignals
import filters
import pyplot
import numpy
from numpy.random import rand

# Just to have fun try this one

def SincBandPassNonLinear(RBTW=0.05):
    """
    Two sinthetic gaussian centered in 15/25 Hz
    example of sinc band pass-filtering
    Relative bandwidth transition desired of 20%
    more info on:
    dspguide.com or wikipedia
    """
    s1 = samplesignals.NonLinear(200, 0.01, 15, 500); # inventa um gauss. centrado em 15Hz
    s2 = samplesignals.NonLinear(200, 0.01, 35, 500); # inventa um gauss. centrado em 35Hz
    s = 0.6*s1+0.4*s2; # monta um com 60% energia de s1 e 40% da energia de s2
    N = numpy.size(s)
    print "Number of samples input %d" % N
    # use the upper limit of frequency to get the number of samples needed to a relative bandwidth transition desired
    Nf = filters.FilterSize(RTbtw=RBTW, Fc=25, dt=0.01);
    fir = filters.SincBandPass(Nf, 0.01, 5, 25); # filtro passa banda caixa frequencias de corte de passagem de 5Hz ah 25Hz
    res = filters.ConvFft(s, fir, 0.01, plot=True); # mostra o resultado da conv.
    # older way
    # res = filters.ConvEnd(s, fir, 0.01, plot=True); # mostra o resultado da conv.
    pyplot.plotfftcompare(s, res, 0.01);
    return res;

def SincTrapezoidalLowPassNoise():
    """
    a experiment with the trapezoidal low pass
    that's perfectly working, filter sampling is perfect
    and also the filtering process with the result
    using NR wrapped around convolution (use FftConv3 its simpler)
    """
    s = samplesignals.Periodic_Noise(200, 0.05);
    N = numpy.size(s)
    print " Input number of samples %d" % numpy.size(s)
    #filter depends on the sample rate and the transition bandwidth we want
    # we want 20% its a reasonable value
    Nf = filters.FilterSize(RTbtw=0.2, Fc=3, dt=0.05);
    print " Filter number of samples %d" % Nf
    filter = filters.SincTrapezoidalLowPass(Nf, 0.05, 0.5, 3);
    #res  = filters.ConvFft(s, filter, 0.05, plot=True);
    res  = filters.ConvFft3(s, filter);
    pyplot.plotfftcompare(s, res, 0.05);

    return res;

def SincTrapezoidalLowPassPureNoise(Dt=1.0, FC=0.3*0.5, Ramp=0.1*0.3*0.5, Signal=rand(100)*10+2, RBTW=0.05):
    """
    testing with noise, rand series values
    you can choose the parameters for the filter
    Dt =sample rate (e.g 1.0)
    Fc = cut-off frequency (e.g 30% of nyquest frequency, just 30% will pass)
    Filter Ramp  = in hertz (e.g 10% of the cutt-of frequency)
    Signal ... (e.g. 100 samples betwen [12, 2])
    RBTW = Relative bandwidth desired of 20%
    Nf = filter Order = Number of samples of the filter is = Nf*2 + 1 (its a Odd number)
    """
    Nf = filters.FilterSize(RBTW, FC, Dt);
    N = numpy.size(Signal)
    print " Input number of samples %d" % N
    #filter smaller than the signal
    filter = filters.SincTrapezoidalLowPass(Nf, Dt, Ramp, FC);
    res = filters.ConvFft(Signal, filter, Dt, plot=True);
    pyplot.plotfftcompare(Signal, res, Dt);

    return res;

def SincBoxNonStationary():
    """
    wrap around convolution from numerical recipes
    including the linear detrend for
    non stationay data
    the fft convolution also makes the things easier for working
    and also performance
    """
    sig = rand(100)+3;
    sig[0:25] = sig[0:25]+1;
    sig[75:100] = sig[75:100]-1;
    filter = filters.SincLowPass(51, 1, 0.1);
    res = filters.ConvFft(sig, filter, 0.1, plot=True);
    pyplot.plotfftcompare(sig, res, 0.1);

def main():
    #SincBandPassNonLinear();
    # more meaningfull and beauty
    SincTrapezoidalLowPassNoise();
    #SincTrapezoidalLowPassPureNoise();
    # good to show solved problem with stationarity using linear detrend
    #SincBoxNonStationary();

if __name__ == '__main__':
    main()

__doc__ = "examples for filters.py"