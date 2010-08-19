#!/usr/bin
# import filters.py

#backend = 'gtk'
#import matplotlib
#matplotlib.use(backend)

import numpy, sys, os
from pylab import arange
import pylab
from pyplot import *
#from numpy import *

import time
import scriptutil # script to manage files etc.. (fonte google)

# Just to have fun try this one
	
def PlayFilter(RBTW=0.05):
	"""
	Two sinthetic gaussian centered in 15/25 Hz 
	example of sinc band pass-filtering
	Relative bandwidth desired of 20% 
	"""
	s1 = sampling2(200, 0.01, 15, 500); # inventa um gauss. centrado em 15Hz
	s2 = sampling2(200, 0.01, 35, 500); # inventa um gauss. centrado em 25Hz	
	s = 0.6*s1+0.4*s2; # monta um com 60% energia de s1 e 40% da energia de s2
	N = numpy.size(s)
	print "Number of samples input %d" % N
	# use the upper limit of frequency to get the number of samples needed to a relative bandwidth transition desired
	Nf = FilterSize(RTbtw=RBTW, Fc=25, dt=0.01);
	fir = boxbanda(Nf, 0.01, 5, 25); # filtro passa banda caixa frequencias de corte de passagem de 5Hz ah 25Hz
	res = resconvFinal(s, fir, 0.01); # mostra o resultado da conv.
	ShowCompleteness(s, res, 0.01);	
	return res;
	
def PlayFilterLowPass():
	"""
	a experiment with the trapezoidal low pass
	that's perfectly working, filter sampling is perfect
	and also the filtering process with the result
	considering the signal as periodic to avoid border effects
	"""
	s = sampling(200, 0.05); 
	N = numpy.size(s)
	print " Input number of samples %d" % numpy.size(s)	
	#filter depends on the sample rate and the transition bandwidth we want
	# we want 20% its a reasonable value
	Nf = FilterSize(RTbtw=0.2, Fc=3, dt=0.05);
	print " Filter number of samples %d" % Nf
	filter = TrapezoidalLowPass(Nf, 0.05, 0.5, 3);
	res  = resconvFinal(s, filter, 0.05);
	ShowCompleteness(s, res, 0.05);	
	
	return res;

def PlayLowPass(Dt=1.0, FC=0.3*0.5, Ramp=0.1*0.3*0.5, Signal=pylab.rand(100)*10+2, RBTW=0.05):
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
	Nf = FilterSize(RBTW, FC, Dt);
	N = numpy.size(Signal)
	print " Input number of samples %d" % N
	#filter smaller than the signal	
	filter = TrapezoidalLowPass(Nf, Dt, Ramp, FC);	
	res = resconvFinal(Signal, filter, Dt);
	ShowCompleteness(Signal, res, Dt);
	
	return res;

def WindowHann(N):
	return numpy.sin(numpy.arange(0, N, 1)*numpy.pi/(N-1));
	

def ShowCompleteness(Sinal, Sinal_Filtered, dt):
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
	plotfftabs(Sinal_Filtered, dt);
	plotfftabs(Sinal, dt);	
	# I dont understand this phase here anyway..
	pylab.subplot(3, 1, 3)
	plotfftphase(Sinal_Filtered, dt);
	plotfftphase(Sinal, dt);
	

# BIG NOTE THE  GUY bellow plot NOT normalized ABS(FFT)
def plotfftabs(signal, dt):
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

def plotfftphase(signal, dt):
	N=numpy.size(signal)
	fund = 1/(N*dt) # freq fundamental, as outras sao multiplas dessa
	f=range(N)
	for k in range(N-1):
		f[k+1]=fund*k #multiplas da fundamental
	z=pylab.fft(signal)
	pylab.plot(f, numpy.angle(z, deg=True))
	pylab.ylabel('Phase (degrees)')
	pylab.xlabel('Frequency (Hz)')

# filtro caixa passa baixa (frequencia corte fcutoff)
# infinite impulse response? e filtro nao causal
# 1) filtro eh truncado no tempo pois nao podemos ter infinitas amostras...
# logo o filtro eh multiplicado por uma funcao caixa no tempo que eh o
# mesmo que convolver com o sinc dessa caixa na frequencia
# que acarreta gibs no modulo e fase?
# 2) o filtro eh amostrado no tempo, ou seja, multiplicado por deltas de
# amplitude 1 e espacados de dt. Que eh o mesmo que convolver na frequencia 
# com deltas espacados de 1/dt e amplitude 1/dt
def box(N, fcutoff, txamos):
	"""
	N is number of samples for the filter
	fcutoff is the F cut-off
	txamos is the sample rate
	Sample the Box Sinc filter simetrically around the zero
	"""			
	
	if(1/(2*txamos) < fcutoff):
		print "Dont be stupid you have to sample the filter at least at the nyquest frequency"
		return
	
	# Amostra simetricamento o operador do filtro em torno do zero
	# nao importa se o numero de pontos for par ou impar, exitem opcoes pros dois casos
	x = y = arange(0, 1, 1)
	
	# caso seja impar amostra incluindo o 0,  amostra perfeitamente simetrico em torno do zero
	#um workaround eh utilizado para evitar divsao por 0 e utilizasse o limite em 0 para setar o valor no 0
	if(N%2!=0): 
		#Nota: arange omite o ultimo valor, por isso o + txamos, para o intervalo ficar [  , ]  e nao [ , )
		x = arange(-txamos*(N-1)/2,(txamos*(N-1)/2)+txamos, txamos) 
		# seta o zero como Nan para evitar, excessao de divisao por 0		
		x[numpy.size(x)/2]=numpy.NaN
		y = txamos*numpy.sin(2*numpy.pi*fcutoff*x)/(numpy.pi*x) # sinc da frequencia de corte, o termo
		# txamos multiplicando serve para garantir a o espc. amplitude em 1, devido a amostragem do filtro no tempo! convolucao com deltas!!
		# inverso da caixa na frequencia
		# txamos multiplicando serve para garantir a o espc. amplitude em 1	
		# set o valor de y no 0, baseado no limite x->0 do operador
		# fazendo o limite (derivando em baixo e em cima) chega-se para Fc => cos(2*pi*Fc*t)*2*Fc
		# com t =0 => 2*Fc, do not forget also the txamos to normalize to 1		
		y[numpy.size(y)/2]= 2*fcutoff*txamos
	else:
		# nao faz muito sentido no mundo real e ninguem usa um filtro nao simetrico mas vai ai
		#amostra  simetricamente em torno do zero sem passar pelo 0, mais facil, nao precisa evitar divisao por zero
		#Nota: arange omite o ultimo valor, por isso o + txamos para o intervalo ficar [  , ]  e nao [ , )
		dispX = txamos*(N.__float__()-1)/2
		x = arange(-dispX,dispX+txamos, txamos) 
		txamos*numpy.sin(2*numpy.pi*fcutoff*x)/(numpy.pi*x); # sinc da frequencia de corte, o termo
		# # txamos multiplicando serve para garantir a o espc. amplitude em 1, devido a amostragem do filtro no tempo! convolucao com deltas!!
	
	
	pylab.figure()
	pylab.subplot(3, 1, 1)
	pylab.plot(x, y)
	pylab.ylabel('Filter sampled')
	print " Filter number of samples %f" % y.__len__()

	pylab.subplot(3, 1, 2)
	plotfftabs(y, txamos);
	pylab.subplot(3, 1, 3)
	plotfftphase(y, txamos);
	
	return y
	
	

def boxbanda(N, txamos, f1, f2):
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
	"""	
	
	if(f2 < f1):
		print "Dont be stupid f2 must be bigger than f1"
		return	
	
	if(1/(2*txamos) < max(f2,f1)):
		print "Dont be stupid you have to sample the filter at least at the nyquest frequency"
		return
	
	Fdelta = (f1+f2)/2 # frequencia central
	Fbox = (f2-f1)/2 # box translated to the central frequency	

	# Amostra simetricamento o operador do filtro em torno do zero
	# nao importa se o numero de pontos for par ou impar, exitem opcoes pros dois casos
	x = y = arange(0, 1, 1)
	
	# caso seja impar amostra incluindo o 0,  amostra perfeitamente simetrico em torno do zero
	#um workaround eh utilizado para evitar divsao por 0 e utilizasse o limite em 0 para setar o valor no 0
	if(N%2!=0): 
		#Nota: arange omite o ultimo valor, por isso o + txamos, para o intervalo ficar [  , ]  e nao [ , )
		x = arange(-txamos*(N-1)/2,(txamos*(N-1)/2)+txamos, txamos) 
		# seta o zero como Nan para evitar, excessao de divisao por 0		
		x[numpy.size(x)/2]=numpy.NaN
		y = 2*txamos*numpy.sin(2*numpy.pi*Fbox*x)*numpy.cos(2*numpy.pi*Fdelta*x)/(numpy.pi*x); # inverso da caixa (soh parte real), o termo
		# txamos multiplicando serve para garantir a o espc. amplitude em 1
		# nao sei o pq do 2? matematicamente		
		# set o valor de y no 0, baseado no limite x->0 do operador
		# fazendo o limite (derivando em baixo e em cima) chega-se para a => cos(2*pi*Fbox*t)*2*Fbox
		# com t =0 => 2*Fbox
		# limite da multplicacao eh a multiplicacao dos limites
		# limite soh do sinc utilizando * 2 anterior
		y[numpy.size(y)/2]= 2*txamos*2*Fbox#2*numpy.pi*Fdelta
		
	else:
		#amostra  simetricamente em torno do zero sem passar pelo 0, mais facil, nao precisa evitar divisao por zero
		#Nota: arange omite o ultimo valor, por isso o + txamos para o intervalo ficar [  , ]  e nao [ , )
		dispX = txamos*(N.__float__()-1)/2
		x = arange(-dispX,dispX+txamos, txamos) 
		y = 2*txamos*numpy.sin(2*numpy.pi*Fbox*x)*numpy.cos(2*numpy.pi*Fdelta*x)/(numpy.pi*x); # inverso da caixa (soh parte real), o termo

	pylab.figure()
	pylab.subplot(3, 1, 1)
	pylab.plot(x, y)
	pylab.ylabel('Filter sampled')
	print " Filter number of samples %f" % y.__len__()
	
	print " Filter Nyquest frequency %f" % (1/(2*txamos))
	print " Filter F1: %f F2: %f" % (f1, f2)
	
	pylab.subplot(3, 1, 2)
	plotfftabs(y, txamos);
	pylab.subplot(3, 1, 3)
	plotfftphase(y, txamos);

	return y


def TrapezoidalLowPass(N, txamos, Ramp, Fc):
	"""
	filtro trapezoidal passa baixa ( rampa = Ramp, Frequencia de corte = Fc )
	infinite impulse response? e filtro nao causal
	1) filtro eh truncado no tempo pois nao podemos ter infinitas amostras...
	logo o filtro eh multiplicado por uma funcao caixa no tempo que eh o
	mesmo que convolver com o sinc dessa caixa na frequencia
	que acarreta gibs no modulo e fase
	2) o filtro eh amostrado no tempo, ou seja, multiplicado por deltas de
	amplitude 1 e espacados de dt. Que eh o mesmo que convolver na frequencia 
	com deltas espacados de 1/dt e amplitude 1/dt
	Txamos = dt  Deve ser  tal que  === 1/2*dt > Fc 
	"""
	a = Ramp/2
	b = Fc - (Ramp/2)
	
	# when sampling the filter operator simetrical to 0
	# the number of samples of the filter is :  N*2+1 where N is the order of the filter
	
	print " Filter Nyquest frequency %f" % (1/(2*txamos))
	print " Filter cut-off frequency Fc: %f Ramp: %f" % (Fc, Ramp)
	
	if(Ramp > Fc):
		print "Dont be stupid Ramp bigger than Fc?"
		return
	
	if(1/(2*txamos) < Fc):
		print "Dont be stupid you have to sample the filter at least at the nyquest frequency"
		return
	
	# this doesnt work, ramp = 0
	if(a==0):
		print "Doesnt make sense a = 0"
		return
	
	# Amostra simetricamento o operador do filtro em torno do zero
	# nao importa se o numero de pontos for par ou impar, exitem opcoes pros dois casos
	x = y = arange(0, 1, 1)
	
	# caso seja impar amostra incluindo o 0,  amostra perfeitamente simetrico em torno do zero
	#um workaround eh utilizado para evitar divsao por 0 e utilizasse o limite em 0 para setar o valor no 0
	if(N%2!=0): 
		#Nota: arange omite o ultimo valor, por isso o + txamos, para o intervalo ficar [  , ]  e nao [ , )
		x = arange(-txamos*(N-1)/2,(txamos*(N-1)/2)+txamos, txamos) 
		# seta o zero como Nan para evitar, excessao de divisao por 0		
		x[numpy.size(x)/2]=numpy.NaN
		y = txamos*numpy.sin(2*numpy.pi*a*x)*numpy.sin(2*numpy.pi*b*x)/(numpy.pi*numpy.pi*x*x*2*a); 
		# inverso da convolucao de
		# duas caixas na frequencia
		# txamos multiplicando serve para garantir a o espc. amplitude em 1	
		# a divisao por 2*a, serve para garantir o spec. amplitude em 1, 
		# pois o resultado da convolucao de duas caixas de amplitude 1 na frequencia (-a, a), (-b, b) com b > a 
		#resulta no trapezio com amplitude maxima igual 2*a		
		# set o valor de y no 0, baseado no limite x->0 do operador
		# fazendo o limite (derivando em baixo e em cima) chega-se para a => cos(2*pi*a*t)*2*a
		# com t =0 => 2*a, ou seja para a e b => 2*a*2*b
		# limite da multplicacao eh a multiplicacao dos limites
		y[numpy.size(y)/2]= 2*a*2*b*txamos/(2*a)
	else:
		#amostra  simetricamente em torno do zero sem passar pelo 0, mais facil, nao precisa evitar divisao por zero
		#Nota: arange omite o ultimo valor, por isso o + txamos para o intervalo ficar [  , ]  e nao [ , )
		dispX = txamos*(N.__float__()-1)/2
		x = arange(-dispX,dispX+txamos, txamos) 
		y = txamos*numpy.sin(2*numpy.pi*a*x)*numpy.sin(2*numpy.pi*b*x)/(numpy.pi*numpy.pi*x*x*2*a);

	pylab.figure()
	pylab.subplot(3, 1, 1)
	pylab.plot(x, y)
	pylab.ylabel('Filter sampled')
	print " Filter number of samples %f" % y.__len__()

	pylab.subplot(3, 1, 2)
	plotfftabs(y, txamos);
	pylab.subplot(3, 1, 3)
	plotfftphase(y, txamos);

	return y
	
	
def sampling1(N=200, dt=0.01): # intervalo de amostragem e amostras 
	#maneira de montar o espectro de frequencia baseado na
	#na expansao em serie de ferier baseado no periodo fundamental
	# dt=0.01; # tx de amostragem
	pylab.figure()
	x0=0 # inicio da amostrag comeca em 0
	x1=dt*N # fim da amostragem
	# k usados para amostrar a funcao
	k=arange(x0,x1,dt) 
	x=k[range(1,N)] # soh N amostras
	fn=1/(2*dt)
	print "frequencia de nyquest = %f" % fn
	print "funcao a amostrar y=2*sin(2*pi*x)+3*cos(2*pi*3*x)+sin(2*pi*3*x)+1.5*cos(2*pi*6*x)"
	y=2*sin(2*pi*x)+3*cos(2*pi*3*x)+sin(2*pi*3*x)+1.5*cos(2*pi*6*x)
	z=fft(y)
	z_real=abs(z)

	fund = 1/(N*dt) # freq fundamental, as outras sao multiplas dessa
	f= range(N-1)
	for k in range(N-2):
		f[k+1]=fund*k
	
	pylab.plot(f)

	#f=0:fn*2/N:2*fn; # o range das frequencias
	pylab.subplot(2,1,1)
	am=range(1,N) # numero de amostras de verdade
	pylab.plot(am,y)
	pylab.xlabel('Amostras')
	pylab.ylabel('Amplitude')
	pylab.subplot(2,1,2)
	pylab.plot(f,z_real);
	pylab.xlabel('frequencia')
	pylab.ylabel('Amplitude')
	
	return (y, am, z_real)

# no dominio continuo, a serie de fourier torna-se transformada de fourier 
# no limite 
#
# frequencia fundamental fp= 1/T com T = (N*dt)
# no maximo N/2*fp


def sampling2(N=200, dt=0.01, fc=50, ct=500): # intervalo de amostragem, amostras, fq. central e constante de decaimento exponencial
	# funcao a amostrar y=exp(c*x^2)
	#na frequencia de nyquist
	# N = 200;
	# dt = 0.01; 50Hz Nyquist
	# fc = 25;
	# ct = 500;
	# espc bonito por exemplo
	pylab.figure()
	i = arange((-N*dt/2),N*dt/2,dt) # divide o mesmo numero de amostras pros dois lados
	x = i[range(1,N)]
	fn=1/(2*dt)
	print "frequencia de nyquest = %f" % fn
	print "funcao a amostrar y=cos(2*pi*fc*x).*exp(-ct*x.*x)"
	y=numpy.cos(2*numpy.pi*fc*x)*numpy.exp(-ct*x*x) #multiplicacao no tempo por um coseno(2*pi*fc*x)
	# igual a convolucao com deltas de frequencia fc
	# igual a deslocar o centro do espectro para a frequencia fc

	z=pylab.fft(y)
	z_real=abs(z)

	pylab.subplot(3,1,1)
	pylab.plot(x,y)
	pylab.ylabel("y=sin(wc*x)*exp(-wc*x^2)")

	pylab.subplot(3,1,2)
	f=arange(0,(2*fn)-(fn*2/N),(fn*2/N)) # o range das frequencias
	pylab.plot(f,z_real)
	pylab.ylabel("Modulo")

	fase=numpy.angle(z, deg=True)
	pylab.subplot(3,1,3)
	pylab.plot(f,fase)
	pylab.ylabel("Fase")
	pylab.xlabel("Frequencia")

	return y

	
def sampling(N, dt): 
	"""
	intervalo de amostragem e amostras
	"""
	#uma maneira de montar a escala do espectro de frequencia baseado
	#na frequencia de nyquist	
	pylab.figure()
	x0=0 # inicio da amostrag comeca em 0
	x1=dt*N # fim da amostragem
	# x usados para amostrar a funcao	
	x=arange(x0,x1,dt) # soh N amostras
	fn=1/(2*dt)
	print "frequencia de nyquest = %f " % fn
	print "funcao a amostrar y=2*sin(2*pi*x)+0.4*cos(2*pi*3*x)+sin(2*pi*3*x)+0.4*cos(2*pi*5*x)"
	y=2*numpy.sin(2*numpy.pi*x)+0.4*numpy.cos(2*numpy.pi*3*x)+numpy.sin(2*numpy.pi*3*x)+0.4*numpy.cos(2*numpy.pi*5*x)
	z=pylab.fft(y)
	z_real=abs(z)
	pylab.subplot(2,1,1)	
	pylab.plot(y)
	pylab.xlabel("Amostras")
	pylab.ylabel("Amplitude")
	pylab.subplot(2,1,2)
	# o range das frequencias
	fund = 1/(N*dt) # freq fundamental, as outras sao multiplas dessa
	f=range(N)
	for k in range(N-1):
		f[k+1]=fund*k #multiplas da fundamental
	pylab.plot(f,z_real)
	pylab.xlabel("frequencia")
	pylab.ylabel("Amplitude")
	return y;


	
def resconvFinal(sinal, filtro, dt):
	"""
	Implemented in C# plug-in for Petrel
	sinal = input signal to filter (time),
	filtro = filter kernel (time) Must be odd number of samples
	dt = sample rate of sinal and filtro
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
	ResConv = pylab.conv(filtro, sinal) # filtra	
	
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
	
	pylab.figure()		
	########################	
	pylab.subplot(5, 2, 1)
	pylab.plot(OrigSignal) #  SINAL A Original TEMPO
	pylab.ylabel("sinal")
	pylab.subplot(5, 2, 2) # SINAL A Original  MODULO Frequencia
	plotfftabs(OrigSignal, dt)
	########################	
	pylab.subplot(5, 2, 3)
	pylab.plot(sinal) #  SINAL A FILTRAR samples added at the begin and at the end, to smooth ending after convolution
	pylab.ylabel("sinal") # a cost that we want
	pylab.subplot(5, 2, 4) #  SINAL A FILTRAR samples added at the begin and at the end Modulo Frequencia
	plotfftabs(sinal, dt)	
	########################
	pylab.subplot(5, 2, 5) 
	pylab.plot(filtro) # FILTRO TEMPO 
	pylab.ylabel("filtro")
	pylab.subplot(5, 2, 6) # FILTRO  MODULO Frequencia
	plotfftabs(filtro, dt)
	########################
	pylab.subplot(5, 2, 7);
	pylab.plot(ResConv) # RESULTADO CONVOLUTION  TEMPO
	pylab.ylabel("Convolution")
	pylab.subplot(5, 2, 8) # RESULTADO CONVOLUTION  MODULO Frequencia
	plotfftabs(ResConv, dt)
	########################
	pylab.subplot(5, 2, 9)
	pylab.plot(ResConvCutted) # RESULTADO CONVOLUTION  CUTTED TEMPO	
	pylab.ylabel('Convolution cutted')	
	pylab.subplot(5, 2, 10) # RESULTADO CONVOLUTION  CUTTED  MODULO Frequencia
	plotfftabs(ResConvCutted, dt)
	###########################

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
	return NextOdd(Fs.__int__());

def NextEven(a):
	if(a%2==1):
		return a+1;
	else:
		return a;

def NextOdd(a):
	if(a%2==0):
		return a+1;
	else:
		return a;

#######################################################	
	################################################
	#########################################
# The code here bellow is result of a fight a great effort to discover
# how to deal with border effects of the convolution 
# and how to clip the result signal in the same size of the input signal
# without changing spectrum or adding undesired effects
# So it's really worthy to remember because I lost A LOT OF TIME WITH THAT

def resconvFFT(sinal, filtro, dt):
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
	sinal = numpy.append(sinal, Zeros(NsRandAdded));

	Ns = numpy.size(sinal) # numero de amostras do sinal a filtrar	
	print " Filtering process Ns %d" % Ns
	
	filtro_extended = filtro
	# to help them to have the same size, N no convolution makes sense otherwise
	if(Nf < Ns):
		# odd - odd = even, so ok, append the begin and the end with zeros
		dif = (Ns - Nf)/2		
		filtro_extended = numpy.append(Zeros(dif), numpy.append(filtro_extended, Zeros(dif)));
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
	#	ResConvCutted = ResConvCutted[range(0,ResConvCutted.__len__()-1)]	
	
	ResConvCutted = ResConv[0:numpy.size(OrigSignal)]
	
	print " Filtering process after cutting process %d" % numpy.size(ResConvCutted)
	print " Filtering process final nyquest %f " % (1/(2*dt))
	
	pylab.figure()	
	
	########################	
	pylab.subplot(6, 2, 1)
	pylab.plot(OrigSignal) #  SINAL A Original TEMPO
	pylab.ylabel("sinal")
	pylab.subplot(6, 2, 2) # SINAL A Original  MODULO Frequencia
	plotfftabs(OrigSignal, dt)
	########################	
	pylab.subplot(6, 2, 3)
	pylab.plot(sinal) #  SINAL A FILTRAR TEMPO EXTENDED
	pylab.ylabel("sinal")
	pylab.subplot(6, 2, 4) # SINAL A FILTRAR  EXTENDED MODULO Frequencia
	plotfftabs(sinal, dt)	
	########################
	pylab.subplot(6, 2, 5) 
	pylab.plot(filtro)  # FILTRO TEMPO 
	pylab.ylabel("filtro")
	pylab.subplot(6, 2, 6) # FILTRO  MODULO Frequencia
	plotfftabs(filtro, dt)		
	########################
	pylab.subplot(6, 2, 7)
	pylab.plot(filtro_extended)  # FILTRO EXTENDED TEMPO
	pylab.ylabel("filtro")
	pylab.subplot(6, 2, 8) # FILTRO ENTENDED MODULO Frequencia
	plotfftabs(filtro_extended, dt)	
	########################
	pylab.subplot(6, 2, 9);
	pylab.plot(ResConv) # RESULTADO CONVOLUTION  TEMPO
	pylab.ylabel("Convolution")
	pylab.subplot(6, 2, 10) # RESULTADO CONVOLUTION  MODULO Frequencia
	plotfftabs(ResConv, dt)
	########################
	pylab.subplot(6, 2, 11)
	pylab.plot(ResConvCutted) # RESULTADO CONVOLUTION  CUTTED TEMPO	
	pylab.ylabel('Convolution cutted')
	# axis([0 length(rscut) min(rscut)-(min(rscut)/50) max(rscut)+(min(rscut)/50) ]);
	pylab.subplot(6, 2, 12) # RESULTADO CONVOLUTION  CUTTED  MODULO Frequencia
	plotfftabs(ResConvCutted, dt)
	###########################

	return  ResConvCutted;

def test():
	"""
	wrap around convolution from numerical recipes
	the next step on implementation, including the linear detrend for
	non stationay data
	the fft convolution also makes the things easier for working
	and also performance
	"""
	sig = pylab.rand(100)+3;
	sig[0:25] = sig[0:25]+1;
	sig[75:100] = sig[75:100]-1;
	filter = box(51, 1, 0.1);
	
	return NRConv(sig, filter, 0.1);


def NRConv(signal, filter, dt, detrend=pylab.detrend_linear):
	"""
	wrap around convolution from numerical recipes
	including the default linear detrend for non stationay data
	the fft convolution also makes the things easier for working (clipping & performance)
	signal : input signal
	filter : must be in simetric around 0
	dt : sample rate
	detrend: detrend function applied before filtering
	"""
	# always detrend the signal with linear detrend if the signal is stationary already it wont be a problem zero mean.
	trend = signal - detrend(signal);
	signal = detrend(signal);
	
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
	fir_hannBeg = numpy.append(fir_hannBeg, numpy.zeros(Ssize-(Fhsize*2+1)));
	
	#put in the wrap around form
	fir_hann = numpy.append(fir_hannBeg, fir_hannEnd);
	
	print "F End %d" % numpy.size(fir_hannEnd)
	print "F Beg %d" % numpy.size(fir_hannBeg)
	print "F %d" % numpy.size(fir_hann)
	print "S %d" % numpy.size(signal)
	
	sigrs = pylab.real(pylab.ifft(pylab.fft(signal)*pylab.fft(fir_hann)));
	
	#python is nice hehhe
	signal = signal[:-(Fhsize+1)] #remove the last part added by the filter avoiding
	sigrs = sigrs[:Ssize] # get the just the first part equivalent to the signal size
	
	pylab.figure();
	pylab.plot(signal+trend);
	pylab.plot(sigrs+trend);
	
	return [signal+trend, sigrs+trend, fir_hann];

def resconv(sinal, filtro, dt):
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
	# that one also didnt work well, the solution found implemented in C# in the plug-in
	# was copying the first and last sample as many times as needed.. that worked perfect
	# needed to be implemented here
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
	
	pylab.figure()		
	########################	
	pylab.subplot(6, 2, 1)
	pylab.plot(OrigSignal) #  SINAL A Original TEMPO
	pylab.ylabel("sinal")
	pylab.subplot(6, 2, 2) # SINAL A Original  MODULO Frequencia
	plotfftabs(OrigSignal, dt)
	########################	
	pylab.subplot(6, 2, 3)
	pylab.plot(sinal) #  SINAL A FILTRAR WRAPPED AROUND TEMPO
	pylab.ylabel("sinal")
	pylab.subplot(6, 2, 4) # SINAL A FILTRAR WRAPPED AROUND MODULO Frequencia
	plotfftabs(sinal, dt)	
	########################
	pylab.subplot(6, 2, 5) 
	pylab.plot(filtro)  # FILTRO TEMPO 
	pylab.ylabel("filtro")
	pylab.subplot(6, 2, 6) # FILTRO  MODULO Frequencia
	plotfftabs(filtro, dt)		
	########################
	pylab.subplot(6, 2, 7)
	pylab.plot(filtro)  # FILTRO EXTENDED TEMPO
	pylab.ylabel("filtro")
	pylab.subplot(6, 2, 8) # FILTRO ENTENDED MODULO Frequencia
	plotfftabs(filtro, dt)	
	########################
	pylab.subplot(6, 2, 9);
	pylab.plot(ResConv) # RESULTADO CONVOLUTION  TEMPO
	pylab.ylabel("Convolution")
	pylab.subplot(6, 2, 10) # RESULTADO CONVOLUTION  MODULO Frequencia
	plotfftabs(ResConv, dt)
	########################
	pylab.subplot(6, 2, 11)
	pylab.plot(ResConvCutted) # RESULTADO CONVOLUTION  CUTTED TEMPO	
	pylab.ylabel('Convolution cutted')
	# axis([0 length(rscut) min(rscut)-(min(rscut)/50) max(rscut)+(min(rscut)/50) ]);
	pylab.subplot(6, 2, 12) # RESULTADO CONVOLUTION  CUTTED  MODULO Frequencia
	plotfftabs(ResConvCutted, dt)
	###########################

	return  ResConvCutted;	
	
	

def ScaledRand(N, min, max):
	return (min+pylab.rand(N)*(max-min))

def Zeros(N):
	return numpy.zeros(N)

def SamplesAroundAddRand(signal, N, min, max):
	"""
	Add N random (between [min, max]) samples at the beginning and at the end
	of the signal array, simetrically
	Signal must have an odd number of samples
	"""
	Ntotal = numpy.size(signal)
	
	#beg = ScaledRand(N, min, max)
	beg = Zeros(N)
	#copy the same to the end, just to not waste time
	end = beg
	
	if(Ntotal%2==0):
		print " Error: to simetrical add samples at beg/end the signal has to have an odd size"
		return 

	return_ = numpy.append(beg, numpy.append(signal, end))	
	
	print " Samples Around: number of rand added %d" % N
	print " Samples Around: final number of samples %d" % numpy.size(return_)
	
	#simetrical samples around 
	return 	return_
	
def SamplesAroundRemove(signal, N):
	"""	
	Remove N samples at the beginning and at the end of the
	signal array, signal must have an odd number of samples	
	"""
	Ntotal = numpy.size(signal)
	
	if(Ntotal%2==0):
		print " Error: to simetrical remove samples at beg/end the signal has to have an odd size"
		return
		
	#Original number of samples, before adding N samples at begin and end
	Noriginal=Ntotal-2*N
	#Noriginal has also an odd number of samples becuase 2*N is always even
	#simetrical samples around Ntotal/2 for Ntotal odd, all Integer operations
	beg = (Ntotal/2)-(Noriginal-1)/2
	end = (Ntotal/2)+(Noriginal-1)/2
	
	#end+1 again because range/arange omits the last one
	# get just the desired samples
	return signal[range(beg,end+1)]

def Normalize(signal):
	"""
	Subtract the mean and 
	divide by the max(signal)-min(signal).
	putting the mean to 0
	and normalizing max and min, betwen 0 and 1
	"""
	tmp = (signal-numpy.mean(signal))
	return tmp/(numpy.max(tmp)-numpy.min(tmp))