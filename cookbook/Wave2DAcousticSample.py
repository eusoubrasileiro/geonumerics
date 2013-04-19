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
import numpy as np

velocity=np.zeros([50,100])
velocity[0:25] = 3000.0
velocity[25:100] = 2500.0
wd2d = Imp2DLuWave(100, 50, 10, 0.001, velocity, 0, 0, 300)
movie= wd2d.Loop()
movie = numpy.save('Imp2DLuWave_100x50x300.npy', movie) # just for preservation save it also
movie = numpy.load('Imp2DLuWave_100x50x300.npy')
WaveAnim.Wave2DAnim(movie, 10, 0.001, velocity, vmin=-0.08, vmax=0.08, filename='Imp2DLuWave_100x50x300')
