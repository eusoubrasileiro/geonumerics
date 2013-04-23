'''
Created on Apr 23, 2013

@author: andre
'''

import numpy
import pylab


def _ConvFft(signal, filterkernel):
    """
    Convolution with fft much faster approach
    works exatcly as convolve(x,y)
    """
    ss = numpy.size(signal);
    fs = numpy.size(filterkernel)
    # padd zeros all until they have the size N+M-1
    signal = numpy.append(signal, numpy.zeros(fs+ss-1-ss));
    filterkernel = numpy.append(filterkernel, numpy.zeros(fs+ss-1-fs));
    signal = pylab.real(pylab.ifft(pylab.fft(signal)*pylab.fft(filterkernel)));
    return signal[:fs+ss-1];

def ConvFft(signal, filterkernel):
    """
    Convolution with fft much faster approach
    works exatcly as convolve(x,y) for y equal a filterkernel response.
    Keeping just the mid part equivalent to the original signal
    so filterkernel must be odd.
    """
    ssor = numpy.size(signal);
    ss = ssor;
    # put one sample 0 more, in case ss is not odd
    # ss+fs-1 end up odd (must be odd to be able to remove the central part convolution)
    #odd+odd -1 = odd
    if(ss%2==0):
        ss = numpy.append(signal, numpy.array([0]));
    ss = numpy.size(signal);
    fs = numpy.size(filterkernel);
    if(fs%2==0):
        print "huahaia get out of here";
        return
    # padd zeros all until they have the size N+M-1
    # and convolve
    signal = _ConvFft(signal, filterkernel)
    beg = (ss+fs-1)/2-(ss-1)/2
    end = (ss+fs-1)/2+(ss-1)/2
    # just the central part of the convolution
    signal = signal[beg-1:end+1]
    # just the original size, avoiding if any zero was added
    return signal[:ssor];


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
    return _NextOdd(Fs.__int__());

def _NextOdd(a):
    if(a%2==0):
        return a+1;
    else:
        return a;

def _Norm(signal):
    """
    Subtract the mean and
    divide by the max(signal)-min(signal).
    putting the mean to 0
    and normalizing max and min, betwen 0 and 1
    """
    tmp = (signal-numpy.mean(signal))
    return tmp/(numpy.max(tmp)-numpy.min(tmp))

# def _LegacyConvFft(signal, filterKernel, dt, detrend=pylab.detrend_linear, plot=False):
#     """
#     wrap around convolution from numerical recipes
#     (wrap around advantage: the only advantage is the the output samples from
#     the convolution are ordered with the desired filtered signal begining in the sample[0])
# 
#     Includes the default linear detrend for non stationay data
#     the fft convolution also makes the things easier for working (clipping & performance)
#     and hanning window over the filterKernel kernel
#     signal : input signal
#     filterKernel : must be in simetric around 0
#     dt : sample rate
#     detrend: detrend function applied before filtering
#     plot = plotting the process or not
#     """
#     # always detrend the signal with linear detrend if the signal is stationary already it wont be a problem zero mean.
#     origsignal = signal
#     trend = signal - detrend(signal);
#     signal = detrend(signal);
#     Sor = numpy.size(signal);
# 
#     if(numpy.size(filterKernel)%2 == 0): #must be odd
#         print "find your way out here huaahhaaa"
#         return
# 
#     Fhsize = (numpy.size(filterKernel)-1)/2; # half size of the filterKernel
#     print "Fhsize %d" % Fhsize
# 
#     signal = numpy.append(signal, numpy.zeros(Fhsize+1)); # easy part just padding zeros, with the half size that is bigger
# 
#     Ssize = numpy.size(signal); # signal size, now signal is bigger, more samples
# 
#     # apply taper to the filterKernel
#     fir_hann = filterKernel*WindowHann(numpy.size(filterKernel));
# 
#     # wrap around form modification
#     # the above creates simetrecally around 0
#     fir_hannEnd = fir_hann[0:Fhsize]; #size ...
#     fir_hannBeg = fir_hann[Fhsize:]; # size ...
#     # make the same size the input signal, considering  the signal size greater than the filterKernel
#     # case the signal is smaller than 2*filterKernel
#     if(Ssize-(Fhsize*2+1) >= 0):
#         fir_hannBeg = numpy.append(fir_hannBeg, numpy.zeros(Ssize-(Fhsize*2+1)));
# 
#     #put in the wrap around form
#     fir_hann = numpy.append(fir_hannBeg, fir_hannEnd);
#     Sfsize = numpy.size(fir_hann);
# 
#     if(Ssize < Sfsize): # the else condition of above
#         #they dont have the same suze
#         # padd the input signal with zeros until they have the same size
#         # to be able to perform the convolution
#         signal = numpy.append(signal, numpy.zeros(Sfsize-Ssize));
# 
#     print "F End %d" % numpy.size(fir_hannEnd)
#     print "F Beg %d" % numpy.size(fir_hannBeg)
#     print "F %d" % numpy.size(fir_hann)
#     print "S %d" % numpy.size(signal)
# 
#     sigrs = pylab.real(pylab.ifft(pylab.fft(signal)*pylab.fft(fir_hann)));
# 
#     #python is nice hehhe
#     # signal = signal[:-(Fhsize+1)] #remove the last part added by the filterKernel avoiding
#     sigrs_ = sigrs[:Sor]; # get the just the first part equivalent to the signal size
#     sigrs_ += trend; #add the trend back
# 
#     if(plot==True):
#         pyplot.plotconvprocess(origsignal, signal, fir_hann, sigrs, sigrs_, dt);
# 
#     return sigrs_;


#######################################################
################################################
#########################################
# # The code here bellow is result of a fight a great effort to discover
# # how to deal with border effects of the convolution
# # and how to clip the result signal in the same size of the input signal
# # without changing spectrum or adding undesired effects
# # So it's really worthy to remember because I lost A LOT OF TIME WITH THAT
# # and all was wrong because the truly problem was about dealing with
# # non stationary data
# 
# def __resconvFFT(sinal, filtro, dt, plot=False):
#     """
#     sinal = input signal to filter (time),
#     filtro = filter kernel (time) Must be odd number of samples
#     dt = sample rate of sinal and filtro
#     calculates the result convolution, and remove the spare samples add by the filter
#     """
# 
#     Ns = numpy.size(sinal) # numero de amostras do sinal a filtrar
#     Nf = numpy.size(filtro) # numero de amostras do filtro
#     OrigSignal = sinal # create a copy of the original signal, just for plotting porpouses
# 
#     if(Nf%2==0):
#         print "The number of samples of the filter must be an odd number"
#         return
# 
#     # to control if a rand sample was added
#     even_samples = 0
# 
#     # para garantir tamanho final identico ao inicial
#     # sem modificacoes na fase, todo sinal de entrada deve ser impar
#     # (if needed) append a rand sample to make it odd. After filtering and  cutting process remove the sample back
#     if(Ns%2==0):
#         even_samples = 1 # we will add one more sample due the signal be even
#         sinal = numpy.append(sinal,  numpy.zeros(1))
#         Ns = numpy.size(sinal)
#         print " New signal size due even number of samples %d" % Ns
# 
#     # Samples to avoid problems with filtering
#     # add the rand samples to avoid boundary effects at the end/beginning
#     NsRandAdded = (Nf-1)/2
#     #sinal = SamplesAroundAddRand(sinal, NsRandAdded, numpy.min(sinal), numpy.max(sinal))
#     sinal = numpy.append(sinal, numpy.zeros(NsRandAdded));
# 
#     Ns = numpy.size(sinal) # numero de amostras do sinal a filtrar
#     print " Filtering process Ns %d" % Ns
# 
#     filtro_extended = filtro
#     # to help them to have the same size, N no convolution makes sense otherwise
#     if(Nf < Ns):
#         # odd - odd = even, so ok, append the begin and the end with zeros
#         dif = (Ns - Nf)/2
#         filtro_extended = numpy.append(numpy.zeros(dif), numpy.append(filtro_extended, numpy.zeros(dif)));
#         Nf = numpy.size(filtro_extended)
#         # arrange in wrap-around mode
#         filtro_extended = numpy.append(filtro_extended[Nf/2:Nf], filtro_extended[0:Nf/2])
# 
#     print " Filtering process Nf %d" % Nf
# 
#     # RESULTADO DA FILTRAGEM
#     ResConv = pylab.ifft(pylab.fft(sinal)*pylab.fft(filtro_extended))
#     ResConv = pylab.real(ResConv )
#     #ResConv = pylab.conv(sinal, filtro_extended) # filtra
#     Nrs = numpy.size(ResConv) # numero de amostras dessa coisa,, Nsinal + Nfiltro - 1 eh impar
# 
#     # RESULTADO DA FILTRAGEM CUTED,
#     # remove o numero total de amostras do filtro
#     # deixa o sinal do tamanho do sinal inicial
#     # Ns  numero de amostras do sinal inicial
#     # Nf numero de amostras do filtro
#     # remove also the rand samples added
# 
#     #beg = (Ns+Nf-1)/2-(Ns-1)/2
#     #end = (Ns+Nf-1)/2+(Ns-1)/2
# 
#     # corta, removendo o numero equivalente as amostras do filtro somente
#     #ResConvCutted = ResConv[range(beg, end+1)]
# 
#     #remove the samples because the Filtering border problems
#     # its not the wisest solution but is what i have
#     #ResConvCutted = SamplesAroundRemove(ResConvCutted, NsRandAdded)
#     #ResConvCutted = ResConvCutted[range(0, ResConvCutted.__len__()-NsRandAdded)]
# 
# 
#     # remove one more at the end
#     #if(even_samples == 1):
#     #    ResConvCutted = ResConvCutted[range(0,ResConvCutted.__len__()-1)]
# 
#     ResConvCutted = ResConv[0:numpy.size(OrigSignal)]
# 
#     print " Filtering process after cutting process %d" % numpy.size(ResConvCutted)
#     print " Filtering process final nyquest %f " % (1/(2*dt))
# 
#     if(plot==True):
#         pyplot.plotconvprocess(OrigSignal, sinal, filtro_extended, ResConv, ResConvCutted)
# 
#     return  ResConvCutted;
# 
# def __resconvOld(sinal, filtro, dt, plot=False):
#     """
#     sinal = input signal to filter (time),
#     filtro = filter kernel (time) Must be odd number of samples
#     dt = sample rate of sinal and filtro
#     calculates the result convolution, and remove the spare samples add by the filter
#     """
# 
#     Ns = numpy.size(sinal) # numero de amostras do sinal a filtrar
#     Nf = numpy.size(filtro) # numero de amostras do filtro
#     OrigSignal = sinal # create a copy of the original signal, just for plotting porpouses
# 
#     if(Nf%2==0):
#         print "The number of samples of the filter must be an odd number"
#         return
# 
#     # to avoid the border effect due filtering
#     # after trying padding with zeros the input signal/ the filter
#     # after put the filter in the wrap-around form
#     # after padding with random values between the min and max values of the input signal
#     # the last solution:
#     # Samples to avoid problems with filtering
#     # add one more Period, do not affect the spectrum.. at the beginning and at the end
#     # creating a new Period in the signal, to workaround the border effect of the convolution
#     # those are signal extensions at the beginning and end,
#     # they contain one period of the signal split in a begin part and end part
#     # that one also didnt work well,
#     Sbegin  = Send = 0;
#     # firt odd number of samples of input signal
#     if(Ns%2==1):
#         #NS is ODDD
#         # copy the b
#         Sbegin = sinal[0 : (Ns-1)/2];
#         Send = sinal[ (Ns-1)/2: Ns];
#     # even case,
#     else:
#         # para garantir tamanho final identico ao inicial
#         # sem modificacoes na fase, todo sinal de entrada deve ser impar
#         # (if needed) append a sample to make it odd. After filtering and  cutting process remove the sample back
#         Sbegin = sinal[0 : (Ns/2)-1];
#         Sbegin = numpy.append(Sbegin, 0);
#         # we also add one more sample due the input signal not being odd, in the begin part that will be putted at the end
#         Send = sinal[ (Ns/2) -1: Ns];
# 
#     # add one more period to the input signal so we can avoid border effects due the convolution
#     sinal = numpy.append(Send, numpy.append(sinal, Sbegin));
# 
#     # RESULTADO DA FILTRAGEM
#     ResConv = numpy.convolve(sinal, filtro) # filtra
# 
#     # RESULTADO DA FILTRAGEM CUTED,
#     # remove o numero total de amostras do filtro
#     # deixa o sinal do tamanho do sinal inicial
#     # Ns  numero de amostras do sinal inicial
#     # Nf numero de amostras do filtro
#     # remove also the rand samples added
#     Ns = numpy.size(sinal) # numero de amostras do sinal a filtrar
#     Nrs = numpy.size(ResConv) # numero de amostras dessa coisa,, Nsinal + Nfiltro - 1 eh impar
# 
#     beg = (Ns+Nf-1)/2-(Ns-1)/2
#     end = (Ns+Nf-1)/2+(Ns-1)/2
# 
#     # corta, removendo o numero equivalente as amostras do filtro somente
#     ResConvCutted = ResConv[range(beg, end+1)]
# 
#     #remove the samples because the Filtering border problems
#     # its the periodic solution
#     print " end %d Nsoriginal %d begin %d" % (numpy.size(Send), numpy.size(OrigSignal), numpy.size(Sbegin))
#     ResConvCutted = ResConvCutted[numpy.size(Send):numpy.size(OrigSignal)+numpy.size(Sbegin)+1]
# 
#     print " Filtering process Ns padded %d" % Ns
#     print " Filtering process after cutting process %d" % numpy.size(ResConvCutted)
#     print " Filtering process final nyquest %f " % (1/(2*dt))
# 
#     if(plot==True):
#         pyplot.plotconvprocess(OrigSignal, sinal, filtro, ResConv, ResConvCutted)
# 
#     return  ResConvCutted;

# def __ConvEndCsharp(sinal, filtro, dt, plot=False):
#     """
#     Implemented in C# plug-in for Petrel v 1.0
#     in the next version changed for a linear detrend to remove non stationarity
# 
#     sinal = input signal to filter (time),
#     filtro = filter kernel (time) Must be odd number of samples
#     dt = sample rate of sinal and filtro
#     to workaround non stationarity, padd the begin and end of the signal
#     with the last sample
#     calculates the result convolution, and remove the spare samples added by the filter
#     """
# 
#     Ns = numpy.size(sinal) # numero de amostras do sinal a filtrar
#     Nf = numpy.size(filtro) # numero de amostras do filtro
#     OrigSignal = numpy.copy(sinal) # create a copy of the original signal, just for plotting porpouses
# 
# 
#     if(Nf%2==0):
#         print "The number of samples of the filter must be an odd number"
#         return
# 
#     # to avoid the border effect due filtering
#     # after trying padding with zeros the input signal/ the filter
#     # after put the filter in the wrap-around form
#     # after padding with random values between the min and max values of the input signal
#     # after add one more Period {
#     # do not affect the spectrum.. at the beginning and at the end
#     # creating a new Period in the signal, to workaround the border effect of the convolution
#     # those are signal extensions at the beginning and end,
#     # they contain one period of the signal split in a begin part and end part
#     #{}
#     #that didnt work well
#     # THE FINAL SOLUTION , the solution found implemented in C# in the plug-in (The phase is not equal the initial in the passing band but who cares??)
#     # was copying the first and last sample as many times as needed.. that worked perfect
#     #the variables bellow they contain begin part (first sample copies) and end part (end sample copies)
#     Sbegin  = Send = 0;
#     # firt odd number of samples of input signal
#     if(Ns%2==1):
#         #NS is ODD, so the begin part is equal the end part
#         # In this case we allways add an EVEN number of samples
#         #at the begin and at the end
#         # Because EVEN+(signal)ODD+EVEN = ODD what we want
#         # the operation bellow just because we need an odd part number of samples
#         NsAdd = Ns/2
#         NsAdd = _NextOdd(NsAdd);
#         Sbegin = numpy.zeros(NsAdd);
#         Sbegin[:] = sinal[0]; # the first sample copies
#         Send = numpy.zeros(NsAdd);
#         Send[:] = sinal[numpy.size(sinal)-1]; # the last sample copies
#     # even case,
#     else:
#         #NS is EVEN, so the begin part is smaller than the end part, because
#         #of one more sample added at the end
#         # Because EVEN+(signal)EVEN+ODD = ODD what we want
#         # added to make the signal odd
#         # the operation bellow just because we need an odd part number of samples
#         NsAdd = Ns/2
#         NsAdd = _NextOdd(NsAdd);
#         Sbegin = numpy.zeros(NsAdd);
#         Sbegin[:] = sinal[0]; # the first sample copies
#         Send = numpy.zeros(NsAdd);
#         Send[:] = sinal[numpy.size(sinal)-1]; # the last sample copies
#         # we also add one more sample due the input signal not being odd, in the end part that will be putted at the end
#         Send = numpy.append(Send, sinal[numpy.size(sinal)-1]);
# 
#     # add one a bunch o samples to the input signal so we can avoid border effects due the convolution
#     # a copied part in the begin and at the end
#     sinal = numpy.append(Sbegin, numpy.append(sinal, Send));
# 
#     # RESULTADO DA FILTRAGEM
#     # convolution is comutative so doesnt matter if it's X (*) Y or Y (*) X
#     ResConv = numpy.convolve(filtro, sinal) # filtra
# 
#     # RESULTADO DA FILTRAGEM CUTED,
#     # remove o numero total de amostras do filtro
#     # deixa o sinal do tamanho do sinal inicial
#     # Ns  numero de amostras do sinal inicial
#     # Nf numero de amostras do filtro
#     # remove also the rand samples added
#     Ns = numpy.size(sinal) # numero de amostras do sinal a filtrar
#     Nrs = numpy.size(ResConv) # numero de amostras dessa coisa,, Nsinal + Nfiltro - 1 eh impar
# 
#     beg = (Ns+Nf-1)/2-(Ns-1)/2
#     end = (Ns+Nf-1)/2+(Ns-1)/2
# 
#     # corta, removendo o numero equivalente as amostras do filtro somente
#     ResConvCutted = ResConv[beg:end];
# 
#     #remove the samples because the Filtering border problems
#     print " Nsoriginal %d begin Added %d end Added %d" % (numpy.size(OrigSignal), numpy.size(Sbegin), numpy.size(Send))
# 
#     beg = numpy.size(Sbegin);
#     end = numpy.size(OrigSignal)+numpy.size(Send);
#     ResConvCutted = ResConvCutted[beg:end];
# 
#     print " Filtering process Ns padded %d" % Ns
#     print " Filtering process after cutting process %d" % numpy.size(ResConvCutted)
#     print " Filtering process final nyquest %f " % (1/(2*dt))
# 
#     if(plot==True):
#         #plot the convolution process
#         pyplot.plotconvprocess(OrigSignal, sinal, filtro, ResConv, ResConvCutted, dt);
# 
#     return  ResConvCutted;
