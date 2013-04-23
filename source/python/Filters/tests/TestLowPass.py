import numpy
from Filters.Sinc import SincLowPass, ConvFft, FilterSize
from numpy.random import rand

"""
 a experiment with the trapezoidal low pass
 that's perfectly working, filterK sampling is perfect
 and also the filtering process with the result
 using NR wrapped around convolution (use FftConv its simpler)
"""

def PeriodicNoise(N=200, dt=0.05):
    """
    periodic signal defined by
    y=2*sin(2*pi*x)+0.4*cos(2*pi*3*x)+sin(2*pi*3*x)+0.4*cos(2*pi*5*x)
    2*2Hz + 0.4*3Hz+ 1*3Hz+0.4*5Hz
    plus random noise
    N = number of samples
    dt = sample rate
    
    """
    #uma maneira de montar a escala do espectro de frequencia baseado
    #na frequencia de nyquist
    x0=0 # inicio da amostrag comeca em 0
    x1=dt*N # fim da amostragem
    # x usados para amostrar a funca
    x=numpy.arange(x0,x1,dt) # soh N amostras
    fn=1/(2*dt)
    print "frequencia de nyquest = %f " % fn

    y=2*numpy.sin(2*numpy.pi*x)+0.4*numpy.cos(2*numpy.pi*3*x)+numpy.sin(2*numpy.pi*3*x)+0.4*numpy.cos(2*numpy.pi*5*x)
    y += -2+4*rand(N);

    return y;


s = PeriodicNoise(200, 0.05);    
print " Input number of samples %d" % numpy.size(s)
#filterK depends on the sample rate and the transition bandwidth we want
# we want 20% its a reasonable value
Nf = FilterSize(RTbtw=0.2, Fc=3, dt=0.05);
print " Filter number of samples %d" % Nf
filterK = SincLowPass(Nf, 3, 0.05);
res  = ConvFft(s, filterK);