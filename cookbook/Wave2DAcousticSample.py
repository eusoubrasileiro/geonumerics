#!/usr/python
"""
runs acoustic wave simulation with 
a simple flat model with 100x50 cells
after save it as an avi animation
"""
import sys
sys.path.append('../source/python');

from WaveFd.Imp2DLuWave import Imp2DLuWave
from WaveFd import WaveAnim
from WaveFd.Wavelet import Triangle
import numpy
velocity=numpy.zeros([50,100])
velocity[0:25] = 3000.0
velocity[25:100] = 2500.0
samplerate=0.001
triangle = Triangle(90.0, samplerate)
wd2d = Imp2DLuWave(100, 50, 10, samplerate, velocity, 0, 0, 700, wavelet=triangle)
movie= wd2d.Simulate()
del wd2d # huge memory usage linear system
movie = numpy.save('Imp2DLuWave_100x50x700.npy', movie) # just for preservation save it also
movie = numpy.load('Imp2DLuWave_100x50x700.npy')
WaveAnim.Wave2DAnim(movie, 10, 0.001, velocity, vmin=-0.002, vmax=0.002, filename='Imp2DLuWave_100x50x700')