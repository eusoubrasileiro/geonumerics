#!/usr/bin
# import filters.py

#backend = 'gtk'
#import matplotlib
#matplotlib.use(backend)

import numpy, sys, os
import pylab
import pyplot
from numpy.random import rand
import time
import scriptutil # script to manage files etc.. (fonte google)


def WindowHann(N):
    """
    Hanning window working for smooth the frequency response
    """
    return numpy.sin(numpy.arange(0, N, 1)*numpy.pi/(N-1));


# filtro caixa passa baixa (frequencia corte fc)
# infinite impulse response? e filtro nao causal
# 1) filtro eh truncado no tempo pois nao podemos ter infinitas amostras...
# logo o filtro eh multiplicado por uma funcao caixa no tempo que eh o
# mesmo que convolver com o sinc dessa caixa na frequencia
# que acarreta gibs no modulo e fase?
# 2) o filtro eh amostrado no tempo, ou seja, multiplicado por deltas de
# amplitude 1 e espacados de dt. Que eh o mesmo que convolver na frequencia
# com deltas espacados de 1/dt e amplitude 1/dt
def SincLowPass(N, fc, dt, plot=False):
    """
    N is number of samples for the filter
    fc is the F cut-off
    dt is the sample rate
    Sample the Box Sinc filter simetrically around the zero
    """

    if(1/(2*dt) < fc):
        print "Dont be stupid you have to sample the filter at least at the nyquest frequency"
        return

    # Amostra simetricamento o operador do filtro em torno do zero
    # nao importa se o numero de pontos for par ou impar, exitem opcoes pros dois casos
    x = y = numpy.arange(0, 1, 1)

    # caso seja impar amostra incluindo o 0,  amostra perfeitamente simetrico em torno do zero
    #um workaround eh utilizado para evitar divsao por 0 e utilizasse o limite em 0 para setar o valor no 0
    if(N%2!=0):
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
    else:
        # nao faz muito sentido no mundo real e ninguem usa um filtro nao simetrico mas vai ai
        #amostra  simetricamente em torno do zero sem passar pelo 0, mais facil, nao precisa evitar divisao por zero
        #Nota: arange omite o ultimo valor, por isso o + dt para o intervalo ficar [  , ]  e nao [ , )
        dispX = dt*(N.__float__()-1)/2
        x = numpy.arange(-dispX,dispX+dt, dt)
        dt*numpy.sin(2*numpy.pi*fc*x)/(numpy.pi*x); # sinc da frequencia de corte, o termo
        # # dt multiplicando serve para garantir a o espc. amplitude em 1, devido a amostragem do filtro no tempo! convolucao com deltas!!

    print " Filter number of samples %f" % y.__len__()

    if(plot==True):
        pyplot.plotfftNabs_phase(y, dt);

    return y



def SincBandPass(N=200, dt=0.01, f1=10, f2=20, plot=False):
    """
    filtro caixa passa banda (frequencias inicias e finais f1 e f2)
    infinite impulse response? e filtro nao causal
    1) filtro eh truncado no tempo pois nao podemos ter infinitas amostras...
    logo o filtro eh multiplicado por uma funcao caixa no tempo que eh o
    mesmo que convolver com o sinc dessa caixa na frequencia
    que acarreta gibs no modulo e fase?
    2) o filtro eh amostrado no tempo, ou seja, multiplicado por deltas de
    amplitude 1 e espacados de dt. Que eh o mesmo que convolver na frequencia
    com deltas espacados de 1/dt e amplitude 1/dt
    To create a band pass filter based on a box, we have to multiply the sinc in time
    for the frequency we want to be the center of our band pass filter.
    Multiplying in time displaces our sinc in frequency, convolving with one delta of the
    that frequency in the frequency domain
    plot = to show or not the response for the filter
    """

    if(f2 < f1):
        print "Dont be stupid f2 must be bigger than f1"
        return

    if(1/(2*dt) < max(f2,f1)):
        print "Dont be stupid you have to sample the filter at least at the nyquest frequency"
        return

    Fdelta = (f1+f2)/2 # frequencia central
    Fbox = (f2-f1)/2 # box translated to the central frequency

    # Amostra simetricamento o operador do filtro em torno do zero
    # nao importa se o numero de pontos for par ou impar, exitem opcoes pros dois casos
    x = y = numpy.arange(0, 1, 1)

    # caso seja impar amostra incluindo o 0,  amostra perfeitamente simetrico em torno do zero
    #um workaround eh utilizado para evitar divsao por 0 e utilizasse o limite em 0 para setar o valor no 0
    if(N%2!=0):
        #Nota: arange omite o ultimo valor, por isso o + dt, para o intervalo ficar [  , ]  e nao [ , )
        x = numpy.arange(-dt*(N-1)/2,(dt*(N-1)/2)+dt, dt)
        # seta o zero como Nan para evitar, excessao de divisao por 0
        x[numpy.size(x)/2]=numpy.NaN
        y = 2*dt*numpy.sin(2*numpy.pi*Fbox*x)*numpy.cos(2*numpy.pi*Fdelta*x)/(numpy.pi*x); # inverso da caixa (soh parte real), o termo
        # dt multiplicando serve para garantir a o espc. amplitude em 1
        # nao sei o pq do 2? matematicamente
        # set o valor de y no 0, baseado no limite x->0 do operador
        # fazendo o limite (derivando em baixo e em cima) chega-se para a => cos(2*pi*Fbox*t)*2*Fbox
        # com t =0 => 2*Fbox
        # limite da multplicacao eh a multiplicacao dos limites
        # limite soh do sinc utilizando * 2 anterior
        y[numpy.size(y)/2]= 2*dt*2*Fbox#2*numpy.pi*Fdelta

    else:
        #amostra  simetricamente em torno do zero sem passar pelo 0, mais facil, nao precisa evitar divisao por zero
        #Nota: arange omite o ultimo valor, por isso o + dt para o intervalo ficar [  , ]  e nao [ , )
        dispX = dt*(N.__float__()-1)/2
        x = numpy.arange(-dispX,dispX+dt, dt)
        y = 2*dt*numpy.sin(2*numpy.pi*Fbox*x)*numpy.cos(2*numpy.pi*Fdelta*x)/(numpy.pi*x); # inverso da caixa (soh parte real), o termo

    print " Filter number of samples %f" % y.__len__()
    print " Filter Nyquest frequency %f" % (1/(2*dt))
    print " Filter F1: %f F2: %f" % (f1, f2)

    if(plot==True):
        pyplot.plotfftNabs_phase(y, dt)

    return y


def SincTrapezoidalLowPass(N, dt, Ramp, Fc, plot=False):
    """
    filtro trapezoidal passa baixa ( rampa = Ramp, Frequencia de corte = Fc )
    infinite impulse response cortada ficando finita (fir) e filtro nao causal
    1) filtro eh truncado no tempo pois nao podemos ter infinitas amostras...
    logo o filtro eh multiplicado por uma funcao caixa no tempo que eh o
    mesmo que convolver com o sinc dessa caixa na frequencia
    que acarreta gibs no modulo e fase
    2) o filtro eh amostrado no tempo, ou seja, multiplicado por deltas de
    amplitude 1 e espacados de dt. Que eh o mesmo que convolver na frequencia
    com deltas espacados de 1/dt e amplitude 1/dt
    dt = dt  Deve ser  tal que  === 1/2*dt > Fc
    """
    a = Ramp/2
    b = Fc - (Ramp/2)

    # when sampling the filter operator simetrical to 0
    # the number of samples of the filter is :  N*2+1 where N is the order of the filter
    # should I change here? the implementation

    print " Filter Nyquest frequency %f" % (1/(2*dt))
    print " Filter cut-off frequency Fc: %f Ramp: %f" % (Fc, Ramp)

    if(Ramp > Fc):
        print "Dont be stupid Ramp bigger than Fc?"
        return

    if(1/(2*dt) < Fc):
        print "Dont be stupid you have to sample the filter at least at the nyquest frequency"
        return

    # this doesnt work, ramp = 0
    if(a==0):
        print "Doesnt make sense a = 0"
        return

    # Amostra simetricamento o operador do filtro em torno do zero
    # nao importa se o numero de pontos for par ou impar, exitem opcoes pros dois casos
    x = y = numpy.arange(0, 1, 1)

    # caso seja impar amostra incluindo o 0,  amostra perfeitamente simetrico em torno do zero
    #um workaround eh utilizado para evitar divsao por 0 e utilizasse o limite em 0 para setar o valor no 0
    if(N%2!=0):
        #Nota: arange omite o ultimo valor, por isso o + dt, para o intervalo ficar [  , ]  e nao [ , )
        x = numpy.arange(-dt*(N-1)/2,(dt*(N-1)/2)+dt, dt)
        # seta o zero como Nan para evitar, excessao de divisao por 0
        x[numpy.size(x)/2]=numpy.NaN
        y = dt*numpy.sin(2*numpy.pi*a*x)*numpy.sin(2*numpy.pi*b*x)/(numpy.pi*numpy.pi*x*x*2*a);
        # inverso da convolucao de
        # duas caixas na frequencia
        # dt multiplicando serve para garantir a o espc. amplitude em 1
        # a divisao por 2*a, serve para garantir o spec. amplitude em 1,
        # pois o resultado da convolucao de duas caixas de amplitude 1 na frequencia (-a, a), (-b, b) com b > a
        #resulta no trapezio com amplitude maxima igual 2*a
        # set o valor de y no 0, baseado no limite x->0 do operador
        # fazendo o limite (derivando em baixo e em cima) chega-se para a => cos(2*pi*a*t)*2*a
        # com t =0 => 2*a, ou seja para a e b => 2*a*2*b
        # limite da multplicacao eh a multiplicacao dos limites
        y[numpy.size(y)/2]= 2*a*2*b*dt/(2*a)
    else:
        #amostra  simetricamente em torno do zero sem passar pelo 0, mais facil, nao precisa evitar divisao por zero
        #Nota: arange omite o ultimo valor, por isso o + dt para o intervalo ficar [  , ]  e nao [ , )
        dispX = dt*(N.__float__()-1)/2
        x = numpy.arange(-dispX,dispX+dt, dt)
        y = dt*numpy.sin(2*numpy.pi*a*x)*numpy.sin(2*numpy.pi*b*x)/(numpy.pi*numpy.pi*x*x*2*a);

    print " Filter number of samples %f" % y.__len__()

    if(plot==True):
        pyplot.plotfftNabs_phase(y, dt)

    return y
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



def ConvFft(signal, filter, dt, detrend=pylab.detrend_linear, plot=False):
    """
    wrap around convolution from numerical recipes
    (wrap around advantage: the only advantage is the the output samples from
    the convolution are ordered with the desired filtered signal begining in the sample[0])

    Includes the default linear detrend for non stationay data
    the fft convolution also makes the things easier for working (clipping & performance)
    and hanning window over the filter kernel
    signal : input signal
    filter : must be in simetric around 0
    dt : sample rate
    detrend: detrend function applied before filtering
    plot = plotting the process or not
    """
    # always detrend the signal with linear detrend if the signal is stationary already it wont be a problem zero mean.
    origsignal = signal
    trend = signal - detrend(signal);
    signal = detrend(signal);
    Sor = numpy.size(signal);

    if(numpy.size(filter)%2 == 0): #must be odd
        print "find your way out here huaahhaaa"
        return

    Fhsize = (numpy.size(filter)-1)/2; # half size of the filter
    print "Fhsize %d" % Fhsize

    signal = numpy.append(signal, numpy.zeros(Fhsize+1)); # easy part just padding zeros, with the half size that is bigger

    Ssize = numpy.size(signal); # signal size, now signal is bigger, more samples

    # apply taper to the filter
    fir_hann = filter*WindowHann(numpy.size(filter));

    # wrap around form modification
    # the above creates simetrecally around 0
    fir_hannEnd = fir_hann[0:Fhsize]; #size ...
    fir_hannBeg = fir_hann[Fhsize:]; # size ...
    # make the same size the input signal, considering  the signal size greater than the filter
    # case the signal is smaller than 2*filter
    if(Ssize-(Fhsize*2+1) >= 0):
        fir_hannBeg = numpy.append(fir_hannBeg, numpy.zeros(Ssize-(Fhsize*2+1)));

    #put in the wrap around form
    fir_hann = numpy.append(fir_hannBeg, fir_hannEnd);
    Sfsize = numpy.size(fir_hann);

    if(Ssize < Sfsize): # the else condition of above
        #they dont have the same suze
        # padd the input signal with zeros until they have the same size
        # to be able to perform the convolution
        signal = numpy.append(signal, numpy.zeros(Sfsize-Ssize));

    print "F End %d" % numpy.size(fir_hannEnd)
    print "F Beg %d" % numpy.size(fir_hannBeg)
    print "F %d" % numpy.size(fir_hann)
    print "S %d" % numpy.size(signal)

    sigrs = pylab.real(pylab.ifft(pylab.fft(signal)*pylab.fft(fir_hann)));

    #python is nice hehhe
    # signal = signal[:-(Fhsize+1)] #remove the last part added by the filter avoiding
    sigrs_ = sigrs[:Sor]; # get the just the first part equivalent to the signal size
    sigrs_ += trend; #add the trend back

    if(plot==True):
        pyplot.plotconvprocess(origsignal, signal, fir_hann, sigrs, sigrs_, dt);

    return sigrs_;

def ConvEnd(sinal, filtro, dt, plot=False):
    """
    Implemented in C# plug-in for Petrel v 1.0
    in the next version changed for a linear detrend to remove non stationarity

    sinal = input signal to filter (time),
    filtro = filter kernel (time) Must be odd number of samples
    dt = sample rate of sinal and filtro
    to workaround non stationarity, padd the begin and end of the signal
    with the last sample
    calculates the result convolution, and remove the spare samples added by the filter
    """

    Ns = numpy.size(sinal) # numero de amostras do sinal a filtrar
    Nf = numpy.size(filtro) # numero de amostras do filtro
    OrigSignal = numpy.copy(sinal) # create a copy of the original signal, just for plotting porpouses


    if(Nf%2==0):
        print "The number of samples of the filter must be an odd number"
        return

    # to avoid the border effect due filtering
    # after trying padding with zeros the input signal/ the filter
    # after put the filter in the wrap-around form
    # after padding with random values between the min and max values of the input signal
    # after add one more Period {
    # do not affect the spectrum.. at the beginning and at the end
    # creating a new Period in the signal, to workaround the border effect of the convolution
    # those are signal extensions at the beginning and end,
    # they contain one period of the signal split in a begin part and end part
    #{}
    #that didnt work well
    # THE FINAL SOLUTION , the solution found implemented in C# in the plug-in (The phase is not equal the initial in the passing band but who cares??)
    # was copying the first and last sample as many times as needed.. that worked perfect
    #the variables bellow they contain begin part (first sample copies) and end part (end sample copies)
    Sbegin  = Send = 0;
    # firt odd number of samples of input signal
    if(Ns%2==1):
        #NS is ODD, so the begin part is equal the end part
        # In this case we allways add an EVEN number of samples
        #at the begin and at the end
        # Because EVEN+(signal)ODD+EVEN = ODD what we want
        # the operation bellow just because we need an odd part number of samples
        NsAdd = Ns/2
        NsAdd = NextEven(NsAdd);
        Sbegin = numpy.zeros(NsAdd);
        Sbegin[:] = sinal[0]; # the first sample copies
        Send = numpy.zeros(NsAdd);
        Send[:] = sinal[numpy.size(sinal)-1]; # the last sample copies
    # even case,
    else:
        #NS is EVEN, so the begin part is smaller than the end part, because
        #of one more sample added at the end
        # Because EVEN+(signal)EVEN+ODD = ODD what we want
        # added to make the signal odd
        # the operation bellow just because we need an odd part number of samples
        NsAdd = Ns/2
        NsAdd = NextEven(NsAdd);
        Sbegin = numpy.zeros(NsAdd);
        Sbegin[:] = sinal[0]; # the first sample copies
        Send = numpy.zeros(NsAdd);
        Send[:] = sinal[numpy.size(sinal)-1]; # the last sample copies
        # we also add one more sample due the input signal not being odd, in the end part that will be putted at the end
        Send = numpy.append(Send, sinal[numpy.size(sinal)-1]);

    # add one a bunch o samples to the input signal so we can avoid border effects due the convolution
    # a copied part in the begin and at the end
    sinal = numpy.append(Sbegin, numpy.append(sinal, Send));

    # RESULTADO DA FILTRAGEM
    # convolution is comutative so doesnt matter if it's X (*) Y or Y (*) X
    ResConv = numpy.convolve(filtro, sinal) # filtra

    # RESULTADO DA FILTRAGEM CUTED,
    # remove o numero total de amostras do filtro
    # deixa o sinal do tamanho do sinal inicial
    # Ns  numero de amostras do sinal inicial
    # Nf numero de amostras do filtro
    # remove also the rand samples added
    Ns = numpy.size(sinal) # numero de amostras do sinal a filtrar
    Nrs = numpy.size(ResConv) # numero de amostras dessa coisa,, Nsinal + Nfiltro - 1 eh impar

    beg = (Ns+Nf-1)/2-(Ns-1)/2
    end = (Ns+Nf-1)/2+(Ns-1)/2

    # corta, removendo o numero equivalente as amostras do filtro somente
    ResConvCutted = ResConv[beg:end];

    #remove the samples because the Filtering border problems
    print " Nsoriginal %d begin Added %d end Added %d" % (numpy.size(OrigSignal), numpy.size(Sbegin), numpy.size(Send))

    beg = numpy.size(Sbegin);
    end = numpy.size(OrigSignal)+numpy.size(Send);
    ResConvCutted = ResConvCutted[beg:end];

    print " Filtering process Ns padded %d" % Ns
    print " Filtering process after cutting process %d" % numpy.size(ResConvCutted)
    print " Filtering process final nyquest %f " % (1/(2*dt))

    if(plot==True):
        #plot the convolution process
        pyplot.plotconvprocess(OrigSignal, sinal, filtro, ResConv, ResConvCutted, dt);

    return  ResConvCutted;

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

def vNorm(signal):
    """
    Subtract the mean and
    divide by the max(signal)-min(signal).
    putting the mean to 0
    and normalizing max and min, betwen 0 and 1
    """
    tmp = (signal-numpy.mean(signal))
    return tmp/(numpy.max(tmp)-numpy.min(tmp))

#######################################################
################################################
#########################################
# The code here bellow is result of a fight a great effort to discover
# how to deal with border effects of the convolution
# and how to clip the result signal in the same size of the input signal
# without changing spectrum or adding undesired effects
# So it's really worthy to remember because I lost A LOT OF TIME WITH THAT
# and all was wrong because the truly problem was about dealing with
# non stationary data

def _resconvFFT(sinal, filtro, dt, plot=False):
    """
    sinal = input signal to filter (time),
    filtro = filter kernel (time) Must be odd number of samples
    dt = sample rate of sinal and filtro
    calculates the result convolution, and remove the spare samples add by the filter
    """

    Ns = numpy.size(sinal) # numero de amostras do sinal a filtrar
    Nf = numpy.size(filtro) # numero de amostras do filtro
    OrigSignal = sinal # create a copy of the original signal, just for plotting porpouses

    if(Nf%2==0):
        print "The number of samples of the filter must be an odd number"
        return

    # to control if a rand sample was added
    even_samples = 0

    # para garantir tamanho final identico ao inicial
    # sem modificacoes na fase, todo sinal de entrada deve ser impar
    # (if needed) append a rand sample to make it odd. After filtering and  cutting process remove the sample back
    if(Ns%2==0):
        even_samples = 1 # we will add one more sample due the signal be even
        sinal = numpy.append(sinal,  Zeros(1))
        Ns = numpy.size(sinal)
        print " New signal size due even number of samples %d" % Ns

    # Samples to avoid problems with filtering
    # add the rand samples to avoid boundary effects at the end/beginning
    NsRandAdded = (Nf-1)/2
    #sinal = SamplesAroundAddRand(sinal, NsRandAdded, numpy.min(sinal), numpy.max(sinal))
    sinal = numpy.append(sinal, numpy.zeros(NsRandAdded));

    Ns = numpy.size(sinal) # numero de amostras do sinal a filtrar
    print " Filtering process Ns %d" % Ns

    filtro_extended = filtro
    # to help them to have the same size, N no convolution makes sense otherwise
    if(Nf < Ns):
        # odd - odd = even, so ok, append the begin and the end with zeros
        dif = (Ns - Nf)/2
        filtro_extended = numpy.append(numpy.zeros(dif), numpy.append(filtro_extended, numpy.zeros(dif)));
        Nf = numpy.size(filtro_extended)
        # arrange in wrap-around mode
        filtro_extended = numpy.append(filtro_extended[Nf/2:Nf], filtro_extended[0:Nf/2])

    print " Filtering process Nf %d" % Nf

    # RESULTADO DA FILTRAGEM
    ResConv = pylab.ifft(pylab.fft(sinal)*pylab.fft(filtro_extended))
    ResConv = pylab.real(ResConv )
    #ResConv = pylab.conv(sinal, filtro_extended) # filtra
    Nrs = numpy.size(ResConv) # numero de amostras dessa coisa,, Nsinal + Nfiltro - 1 eh impar

    # RESULTADO DA FILTRAGEM CUTED,
    # remove o numero total de amostras do filtro
    # deixa o sinal do tamanho do sinal inicial
    # Ns  numero de amostras do sinal inicial
    # Nf numero de amostras do filtro
    # remove also the rand samples added

    #beg = (Ns+Nf-1)/2-(Ns-1)/2
    #end = (Ns+Nf-1)/2+(Ns-1)/2

    # corta, removendo o numero equivalente as amostras do filtro somente
    #ResConvCutted = ResConv[range(beg, end+1)]

    #remove the samples because the Filtering border problems
    # its not the wisest solution but is what i have
    #ResConvCutted = SamplesAroundRemove(ResConvCutted, NsRandAdded)
    #ResConvCutted = ResConvCutted[range(0, ResConvCutted.__len__()-NsRandAdded)]


    # remove one more at the end
    #if(even_samples == 1):
    #    ResConvCutted = ResConvCutted[range(0,ResConvCutted.__len__()-1)]

    ResConvCutted = ResConv[0:numpy.size(OrigSignal)]

    print " Filtering process after cutting process %d" % numpy.size(ResConvCutted)
    print " Filtering process final nyquest %f " % (1/(2*dt))

    if(plot==True):
        pyplot.plotconvprocess(OrigSignal, sinal, filtro_extended, ResConv, ResConvCutted)

    return  ResConvCutted;

def _resconvOld(sinal, filtro, dt, plot=False):
    """
    sinal = input signal to filter (time),
    filtro = filter kernel (time) Must be odd number of samples
    dt = sample rate of sinal and filtro
    calculates the result convolution, and remove the spare samples add by the filter
    """

    Ns = numpy.size(sinal) # numero de amostras do sinal a filtrar
    Nf = numpy.size(filtro) # numero de amostras do filtro
    OrigSignal = sinal # create a copy of the original signal, just for plotting porpouses

    if(Nf%2==0):
        print "The number of samples of the filter must be an odd number"
        return

    # to avoid the border effect due filtering
    # after trying padding with zeros the input signal/ the filter
    # after put the filter in the wrap-around form
    # after padding with random values between the min and max values of the input signal
    # the last solution:
    # Samples to avoid problems with filtering
    # add one more Period, do not affect the spectrum.. at the beginning and at the end
    # creating a new Period in the signal, to workaround the border effect of the convolution
    # those are signal extensions at the beginning and end,
    # they contain one period of the signal split in a begin part and end part
    # that one also didnt work well,
    Sbegin  = Send = 0;
    # firt odd number of samples of input signal
    if(Ns%2==1):
        #NS is ODDD
        # copy the b
        Sbegin = sinal[0 : (Ns-1)/2];
        Send = sinal[ (Ns-1)/2: Ns];
    # even case,
    else:
        # para garantir tamanho final identico ao inicial
        # sem modificacoes na fase, todo sinal de entrada deve ser impar
        # (if needed) append a sample to make it odd. After filtering and  cutting process remove the sample back
        Sbegin = sinal[0 : (Ns/2)-1];
        Sbegin = numpy.append(Sbegin, 0);
        # we also add one more sample due the input signal not being odd, in the begin part that will be putted at the end
        Send = sinal[ (Ns/2) -1: Ns];

    # add one more period to the input signal so we can avoid border effects due the convolution
    sinal = numpy.append(Send, numpy.append(sinal, Sbegin));

    # RESULTADO DA FILTRAGEM
    ResConv = pylab.conv(sinal, filtro) # filtra

    # RESULTADO DA FILTRAGEM CUTED,
    # remove o numero total de amostras do filtro
    # deixa o sinal do tamanho do sinal inicial
    # Ns  numero de amostras do sinal inicial
    # Nf numero de amostras do filtro
    # remove also the rand samples added
    Ns = numpy.size(sinal) # numero de amostras do sinal a filtrar
    Nrs = numpy.size(ResConv) # numero de amostras dessa coisa,, Nsinal + Nfiltro - 1 eh impar

    beg = (Ns+Nf-1)/2-(Ns-1)/2
    end = (Ns+Nf-1)/2+(Ns-1)/2

    # corta, removendo o numero equivalente as amostras do filtro somente
    ResConvCutted = ResConv[range(beg, end+1)]

    #remove the samples because the Filtering border problems
    # its the periodic solution
    print " end %d Nsoriginal %d begin %d" % (numpy.size(Send), numpy.size(OrigSignal), numpy.size(Sbegin))
    ResConvCutted = ResConvCutted[numpy.size(Send):numpy.size(OrigSignal)+numpy.size(Sbegin)+1]

    print " Filtering process Ns padded %d" % Ns
    print " Filtering process after cutting process %d" % numpy.size(ResConvCutted)
    print " Filtering process final nyquest %f " % (1/(2*dt))

    if(plot==True):
        pyplot.plotconvprocess(OrigSignal, sinal, filtro, ResConv, ResConvCutted)

    return  ResConvCutted;