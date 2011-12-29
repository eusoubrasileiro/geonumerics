import numpy, sys, os
import pylab
from numpy.random import rand
import time

def WindowHann(N):
    """
    Hanning window working for smooth the frequency response
    """
    return numpy.sin(numpy.arange(0, N, 1)*numpy.pi/(N-1));

def ConvFft2(signal, filter):
    """
    Convolution with fft much faster approach
    works exatcly as convolve(x,y)
    """
    ss = numpy.size(signal);
    fs = numpy.size(filter)
    # padd zeros all until they have the size N+M-1
    signal = numpy.append(signal, numpy.zeros(fs+ss-1-ss));
    filter = numpy.append(filter, numpy.zeros(fs+ss-1-fs));
    signal = pylab.real(pylab.ifft(pylab.fft(signal)*pylab.fft(filter)));
    return signal[:fs+ss-1];

def ConvFft3(signal, filter):
    """
    Convolution with fft much faster approach
    works exatcly as convolve(x,y) for y equal a filter response.
    Conserving just the mid part equivalent to the original signal
    so filter must be odd.
    """
    ssor = numpy.size(signal);
    ss = ssor;
    # put one sample 0 more, in case ss is not odd
    # ss+fs-1 end up odd (must be odd to be able to remove the central part convolution)
    #odd+odd -1 = odd
    if(ss%2==0):
        ss = numpy.append(signal, numpy.array([0]));
    ss = numpy.size(signal);
    fs = numpy.size(filter);
    if(fs%2==0):
        print "huahaia get out of here";
        return
    # padd zeros all until they have the size N+M-1
    # and convolve
    signal = ConvFft2(signal, filter)
    beg = (ss+fs-1)/2-(ss-1)/2
    end = (ss+fs-1)/2+(ss-1)/2
    # just the central part of the convolution
    signal = signal[beg-1:end+1]
    # just the original size, avoiding if any zero was added
    return signal[:ssor];

# filtro caixa passa baixa (frequencia corte fc)
# infinite impulse response? e filtro nao causal
# 1) filtro eh truncado no tempo pois nao podemos ter infinitas amostras...
# logo o filtro eh multiplicado por uma funcao caixa no tempo que eh o
# mesmo que convolver com o sinc dessa caixa na frequencia
# que acarreta gibs no modulo e fase?
# 2) o filtro eh amostrado no tempo, ou seja, multiplicado por deltas de
# amplitude 1 e espacados de dt. Que eh o mesmo que convolver na frequencia
# com deltas espacados de 1/dt e amplitude 1/dt
def SincLowPass(N, fc, dt):
    """
    N is number of samples for the filter always odd
    fc is the F cut-off
    dt is the sample rate
    Sample the Box Sinc filter simetrically around the zero
    """
    if(N%2==0):
        return

    if(1/(2*dt) < fc):
        print "Dont be stupid you have to sample the filter at least at the nyquest frequency"
        return

    # Amostra simetricamento o operador do filtro em torno do zero
    # nao importa se o numero de pontos for par ou impar, exitem opcoes pros dois casos
    x = y = numpy.arange(0, 1, 1)

    #Nota: arange omite o ultimo valor, por isso o + dt, para o intervalo ficar [  , ]  e nao [ , )
    x = numpy.arange(-dt*(N-1)/2,(dt*(N-1)/2)+dt, dt)
    # seta o zero como Nan para evitar, excessao de divisao por 0
    x[numpy.size(x)/2]=numpy.NaN
    y = dt*numpy.sin(2*numpy.pi*fc*x)/(numpy.pi*x) # sinc da frequencia de corte, o termo
    # dt multiplicando serve para garantir a o espc. amplitude em 1, devido a amostragem do filtro no tempo! convolucao com deltas!!
    # inverso da caixa na frequencia
    # dt multiplicando serve para garantir a o espc. amplitude em 1
    # set o valor de y no 0, baseado no limite x->0 do operador
    # fazendo o limite (derivando em baixo e em cima) chega-se para Fc => cos(2*pi*Fc*t)*2*Fc
    # com t =0 => 2*Fc, do not forget also the dt to normalize to 1
    y[numpy.size(y)/2]= 2*fc*dt   

    print " Filter number of samples %f" % y.__len__()

    return y



def FilterSize(RTbtw=0.2, Fc=125, dt=0.001):
    """
    Gives the number of samples (odd number) needed to sample our filter operator
    based on:
    1) the relative transition bandwidth RTbtw WE WANT
    (Transition band size divide by cut-off frequency)
    2) Filter cut-off frequency
    3) Sample rate that's gonna be used to sample the filter
    This aproximation is valid for Sinc Filters.
    a RTbtw = 20% is a reasonable value for a non distorded frequency response
    """
    Fs = 2/(RTbtw*Fc*dt);
    return nextodd(Fs.__int__());

def nextodd(a):
    if(a%2==0):
        return a+1;
    else:
        return a;

def Periodic_Noise(N=200, dt=0.05):
    """
    periodic signal defined by
    y=2*sin(2*pi*x)+0.4*cos(2*pi*3*x)+sin(2*pi*3*x)+0.4*cos(2*pi*5*x)
    2*2Hz + 0.4*3Hz+ 1*3Hz+0.4*5Hz
    plus random noise
    N = number of samples
    dt = sample rate
    plot = show espectrum after
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

def SincLowPassNoise():
    """
    a experiment with the trapezoidal low pass
    that's perfectly working, filter sampling is perfect
    and also the filtering process with the result
    using NR wrapped around convolution (use FftConv3 its simpler)
    """
    s = Periodic_Noise(200, 0.05);
    N = numpy.size(s)
    print " Input number of samples %d" % numpy.size(s)
    #filter depends on the sample rate and the transition bandwidth we want
    # we want 20% its a reasonable value
    Nf = FilterSize(RTbtw=0.2, Fc=3, dt=0.05);
    print " Filter number of samples %d" % Nf
    filter = SincLowPass(Nf, 3, 0.05);
    res  = ConvFft3(s, filter);

    return [res, s, filter];