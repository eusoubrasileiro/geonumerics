#!/usr/python
"""
runs acoustic wave simulation with 
a simple flat model with 100x50 cells
after save it as an avi animation
"""
import sys
sys.path.append('../source/python');

from WaveFd.Imp2DLuWave import Imp2DLuSparseWave
from WaveFd import WaveAnim
from WaveFd.Wavelet import RickerSource
import numpy
velocity=numpy.zeros([50,100])
velocity[0:25] = 3000.0
velocity[25:100] = 2500.0
samplerate=0.001
ricker = RickerSource(90.0, 1.0, 1./90.0)
wd2d = Imp2DLuSparseWave(100, 50, 10., samplerate, velocity, 0, 0, 700, wavelet=ricker)
movie= wd2d.Simulate()
del wd2d # huge memory usage linear system
movie = numpy.save('Imp2DLuWave_100x50x700.npy', movie) # just for preservation save it also
movie = numpy.load('Imp2DLuWave_100x50x700.npy')
WaveAnim.Wave2DAnim(movie, 10, samplerate, velocity, filename='Imp2DLuWave_100x50x700')