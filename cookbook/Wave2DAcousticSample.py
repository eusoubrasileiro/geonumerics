#!/usr/python
"""
runs acoustic wave simulation with 
a simple flat model with 100x100 cells
after save it as an avi animation
"""
sys.path.append('../source/python');

from WaveFd.Imp2DLuWave import Imp2DLuWave
from WaveFd import WaveAnim, Wavelet
import numpy

velocity=zeros([100,100])
velocity[0:50] = 3000.0
velocity[50:100] = 4000.0
wavelet = 10*Wavelet.Triangle(Fc=125,Dt=0.001) # create a custom wavelet
wd2d = Imp2DLuWave(100, 100, 10, 0.001, velocity, 0, 0, 500, wavelet=wavelet)
movie= wd2d.Loop()
movie = numpy.save('Imp2DLuWave_100x100x500_125_Fc.npy', movie) # just for preservation save it also
movie = numpy.load('Imp2DLuWave_100x100x500_125_Fc.npy')
WaveAnim.Wave2DAnim(movie, 10, 0.001, velocity, vmin=-0.05, vmax=0.05, filename='Imp2DLuWave_100x100x500_125_Fc')
