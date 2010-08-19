#!/usr/bin
#script python pra plotar 
#para usar vc tem que fazer
# import plotfile
#e depois
# plotfile.plotfile('sq.txt', 'b^-')

import numpy, sys, os
import time
import scriptutil # script to manage files etc.. (fonte google)
import pylab

__doc__ = "These are fuctions to helping on plotting/loading you can also use the load, save from pylab"

def showfile(filename, opt=None, lwidth=1.5):
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
	

def plotfile(filename, opt=None, lwidth=1.5): 
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

def plotPowersSpectrumdB(y, dt=0.1):
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
	
def plotfftNabs(y, dt=0.1): # plot the grah of fft (amplitude) of y with command plot
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
	
def plotfftphaseR(y, dt=0.1): # plot the grah of fft (fase)  of y with command plot
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
	