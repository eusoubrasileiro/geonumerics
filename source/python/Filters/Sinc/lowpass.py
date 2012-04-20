import numpy
import pylab

"""
To play use SincLowPass()
to create the filter kernel 
and ConvFft() to apply it

Eg.

lowpass.ConvFft(Sd, lowpass.SincLowPass(1001, 0.5, 0.1) )

filters the Sd array using a LowPass filter of 1001 points
with a cut-off frequency of 0.5 Hz and sample rate of 0.1 s

You can also use a hanning taper window 
before applying the filter

Eg.
fkernel = lowpass.SincLowPass(1001, 0.5, 0.1) 
fkernel_tapered = Fkernel*WindowHann(np.size(fkernel))
...

"""

def SincLowPass(N, fc, dt):
    """
    N is number of samples for the filter always odd
    fc is the F cut-off
    dt is the sample rate
    Sample the Box Sinc filter simetrically around the zero
     filtro caixa passa baixa (frequencia corte fc)
    infinite impulse response? e filtro nao causal
    1) filtro eh truncado no tempo pois nao podemos ter infinitas amostras...
    logo o filtro eh multiplicado por uma funcao caixa no tempo que eh o
    mesmo que convolver com o sinc dessa caixa na frequencia
    que acarreta gibs no modulo e fase?
    2) o filtro eh amostrado no tempo, ou seja, multiplicado por deltas de
    amplitude 1 e espacados de dt. Que eh o mesmo que convolver na frequencia
    com deltas espacados de 1/dt e amplitude 1/dt
    """
    if(N%2==0):
        return

    if(1/(2*dt) < fc):
        print "Dont be stupid you have to sample the filter at least at the nyquest frequency"
        return

    # Amostra simetricamento o operador do filtro em torno do zero
    # nao importa se o numero de pontos for par ou impar, exitem opcoes pros dois casos    

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


def WindowHann(N):
    """
    Hanning window working for smooth the frequency response
    """
    return numpy.sin(numpy.arange(0, N, 1)*numpy.pi/(N-1));

def _NextOdd(a):
    if(a%2==0):
        return a+1;
    else:
        return a;

def _FilterSize(RTbtw=0.2, Fc=125, dt=0.001):
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
    return _NextOdd(Fs.__int__());

def _ConvFft(signal, FilterKernel): 
    """
    Convolution with fft much faster approach
    works exactly as convolve(x,y)
    """
    ss = numpy.size(signal);
    fs = numpy.size(FilterKernel)
    # padd zeros all until they have the size N+M-1
    signal = numpy.append(signal, numpy.zeros(fs+ss-1-ss));
    FilterKernel = numpy.append(FilterKernel, numpy.zeros(fs+ss-1-fs));
    signal = pylab.real(pylab.ifft(pylab.fft(signal)*pylab.fft(FilterKernel)));
    return signal[:fs+ss-1];

def ConvFft(signal, FilterKernel):
    """
    Convolution with fft much faster approach
    works exatcly as convolve(x,y) for y equal a FilterKernel response.
    Conserving just the mid part equivalent to the original signal
    so FilterKernel must be odd.
    """
    ssor = numpy.size(signal);
    ss = ssor;
    # put one sample 0 more, in case ss is not odd
    # ss+fs-1 end up odd (must be odd to be able to remove the central part convolution)
    #odd+odd -1 = odd
    if(ss%2==0):
        ss = numpy.append(signal, numpy.array([0]));
    ss = numpy.size(signal);
    fs = numpy.size(FilterKernel);
    if(fs%2==0):
        print "huahaia get out of here";
        return
    # padd zeros all until they have the size N+M-1
    # and convolve
    signal = _ConvFft(signal, FilterKernel)
    beg = (ss+fs-1)/2-(ss-1)/2
    end = (ss+fs-1)/2+(ss-1)/2
    # just the central part of the convolution
    signal = signal[beg-1:end+1]
    # just the original size, avoiding if any zero was added
    return signal[:ssor];







