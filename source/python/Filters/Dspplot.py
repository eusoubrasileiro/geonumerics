#!/usr/bin
#script python pra plotar
#para usar vc tem que fazer
# import plotfile
#e depois
# pyplot.plotfile('sq.txt', 'b^-')

import numpy
import pylab

"""
These are fuctions to helping on plotting/loading 
you can also use the load, save from pylab
"""

def ShowFile(filename, opt=None, lwidth=1.5):
    """
    Plot the archive filename
    and after use the command show()
    """
    f = pylab.load(filename)

    if numpy.size(numpy.shape(f)) == 1: # uma dimensao o arquivo
        if opt:
            pylab.plot(range(0, len(f), 1), f, opt, linewidth=lwidth)
        else :
            pylab.plot(range(0, len(f), 1), f, linewidth=lwidth)

    if numpy.size(numpy.shape(f)) == 2: # duas dimensoes o arquivo
        if opt:
            pylab.plot(f[:,0],f[:,1], opt, linewidth=lwidth)
        else :
            pylab.plot(f[:,0],f[:,1], linewidth=lwidth)

    pylab.show();

    return f


def PlotFile(filename, opt=None, lwidth=1.5):
    """
    Plot the archive filename
    without the command show()
    """
    f = pylab.load(filename)

    if numpy.size(numpy.shape(f)) == 1: # uma dimensao o arquivo
        if opt :
            pylab.plot(range(0, len(f), 1), f, opt, linewidth=lwidth)
        else :
            pylab.plot(range(0, len(f), 1), f, linewidth=lwidth)

    if numpy.size(numpy.shape(f)) == 2: # duas dimensoes o arquivo
        if opt :
            pylab.plot(f[:,0],f[:,1], opt, linewidth=lwidth)
        else :
            pylab.plot(f[:,0],f[:,1], linewidth=lwidth)

    return f

def PlotPowersSpectrumdB(y, dt=0.1):
    """
    Try of power spectral density estimation
    not windowed no Welch yet.
    plot also just half spectrum (just positive frequencies),
    and from frequencys 0:N/2
    the mean (DC component) is removed before fft
    """
    # plot the grah of fft (amplitude) of y with command plot
    N = numpy.size(y)

    #remove the DC component just the mean (could use a detrend function)
    y = y - numpy.mean(y)

    fund = 1/(N*dt) # freq fundamental, as outras sao multiplas dessa
    f=range(N)
    for k in range(N-1):
        f[k+1]=fund*k #multiplas da fundamental

    yn = pylab.fft(y)
    power = yn*numpy.conj(yn)/(numpy.sqrt(numpy.size(yn)))
    power = numpy.abs(power)
    # create decibel scale Lb = log10(P1/P0) where P0 is the maximum value
    P0 = max(power)
    power = numpy.log10(power/P0);

    #ignore Dc and the negative frequencies
    pylab.plot(f[0:numpy.size(f)/2], power[0:numpy.size(f)/2])

def PlotFftNAbs(y, dt=0.1): # plot the grah of fft (amplitude) of y with command plot
    """
    Normalized version plot of FFT
    equal abs(FFT) * 2/y.length
    plot also just half spectrum (just positive frequencies),
    and from frequencys 1:N/2
    the mean (DC component) is removed before fft
    """
    N = numpy.size(y)

    #remove the DC component just the mean (could use a detrend function)
    y = y - numpy.mean(y)

    fund = 1/(N*dt) # freq fundamental, as outras sao multiplas dessa
    f=range(N)
    for k in range(N-1):
        f[k+1]=fund*k #multiplas da fundamental

    yn = pylab.fft(y)
    # perfeito. se vpce colocar um seno de 30 e outro de 60, com aplitude 1 e 2 no espectro
    # de amplitude voce vai encontrar exatamente 1 e 2
    absyn=2*numpy.abs(yn)/numpy.size(yn);
    pylab.plot(f[0:numpy.size(f)/2], absyn[0:numpy.size(f)/2])
    return [f, absyn]

# BIG NOTE THE  GUY bellow plot NOT normalized ABS(FFT)
def PlotFftAbs(signal, dt):
    N=numpy.size(signal)
    fund = 1/(N*dt) # freq fundamental, as outras sao multiplas dessa
    f=range(N)
    for k in range(N-1):
        f[k+1]=fund*k #multiplas da fundamental
    z=pylab.fft(signal)
    z_real=abs(z)
    pylab.plot(f, z_real)
    pylab.ylabel('Amplitude')
    pylab.xlabel('Frequency (Hz)')

def PlotFftPhase(signal, dt):
    N=numpy.size(signal)
    fund = 1/(N*dt) # freq fundamental, as outras sao multiplas dessa
    f=range(N)
    for k in range(N-1):
        f[k+1]=fund*k #multiplas da fundamental
    z=pylab.fft(signal)
    pylab.plot(f, numpy.angle(z, deg=True))
    pylab.ylabel('Phase (degrees)')
    pylab.xlabel('Frequency (Hz)')

def PlotFftRPhase(y, dt=0.1): # plot the grah of fft (fase)  of y with command plot
    """
    plot also just half spectrum (just positive frequencies),
    and from frequencys 1:N/2
    """
    N = numpy.size(y)

    fund = 1/(N*dt) # freq fundamental, as outras sao multiplas dessa
    f=range(N)
    for k in range(N-1):
        f[k+1]=fund*k #multiplas da fundamental

    yn = pylab.fft(y)
    angyn = numpy.angle(yn, deg=True)

    pylab.plot(f[0:numpy.size(f)/2], angyn[0:numpy.size(f)/2])


def PlotFftNAbsPhase(y, dt, ploty=False):
    """
    plot the fft abs spectrum and phase
    if ploty= true also plot the input signal
    """
    if(ploty==True):
        pylab.subplot(3,1,1);
        pylab.plot(y);
        pylab.subplot(3,1,2);
        PlotFftNAbs(y, dt);
        pylab.subplot(3,1,3);
        PlotFftRPhase(y, dt);
    else:
        pylab.subplot(2,1,1);
        PlotFftNAbs(y, dt);
        pylab.subplot(2,1,2);
        PlotFftRPhase(y, dt);


def PlotFftCompare(Sinal, Sinal_Filtered, dt):
    """
    no idea of how to implement this... how to show completeness..?
    just plotting comparisons between a filtered signal and its original version
    """
    pylab.figure();
    pylab.subplot(3, 1, 1) # all togueter 3 in the same graph
    pylab.plot(Sinal_Filtered);
    pylab.plot(Sinal);
    pylab.xlabel('Time')
    pylab.subplot(3, 1, 2)
    PlotFftAbs(Sinal_Filtered, dt);
    PlotFftAbs(Sinal, dt);
    # I dont understand this phase here anyway..
    pylab.subplot(3, 1, 3)
    PlotFftPhase(Sinal_Filtered, dt);
    PlotFftPhase(Sinal, dt);

def PlotConvProcess(origsignal, paddedsignal, filterKernel, resconv, resconvcutted, dt):
    """
    Plot the convolution process, time & frequency
    parameters:
    original signal
    padded signal before effective convolution
    filterKernel kernel
    signal after convolution
    signal after convolution samples clipped to maintain the original size
    dt = sample rate of all guys
    """
    pylab.figure()
    ########################
    pylab.subplot(5, 2, 1)
    pylab.plot(origsignal) #  SINAL A Original TEMPO
    pylab.ylabel("signal")
    pylab.subplot(5, 2, 2) # SINAL A Original  MODULO Frequencia
    PlotFftAbs(origsignal, dt)
    ########################
    pylab.subplot(5, 2, 3)
    pylab.plot(paddedsignal) #  SINAL A FILTRAR samples added at the begin and at the end, to smooth ending after convolution
    pylab.ylabel("signal") # a cost that we want
    pylab.subplot(5, 2, 4) #  SINAL A FILTRAR samples added at the begin and at the end Modulo Frequencia
    PlotFftAbs(paddedsignal, dt)
    ########################
    pylab.subplot(5, 2, 5)
    pylab.plot(filterKernel) # FILTRO TEMPO
    pylab.ylabel("filterKernel")
    pylab.subplot(5, 2, 6) # FILTRO  MODULO Frequencia
    PlotFftAbs(filterKernel, dt)
    ########################
    pylab.subplot(5, 2, 7);
    pylab.plot(resconv) # RESULTADO CONVOLUTION  TEMPO
    pylab.ylabel("Convolution")
    pylab.subplot(5, 2, 8) # RESULTADO CONVOLUTION  MODULO Frequencia
    PlotFftAbs(resconv, dt)
    ########################
    pylab.subplot(5, 2, 9)
    pylab.plot(resconvcutted) # RESULTADO CONVOLUTION  CUTTED TEMPO
    pylab.ylabel('Convolution cutted')
    pylab.subplot(5, 2, 10) # RESULTADO CONVOLUTION  CUTTED  MODULO Frequencia
    PlotFftAbs(resconvcutted, dt)
    ###########################