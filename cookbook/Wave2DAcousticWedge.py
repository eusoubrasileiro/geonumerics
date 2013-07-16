#!/usr/python
"""
runs acoustic wave simulation with 
a simple wedge model with 400x160 cells.
4000 bellow 2000 above velocities 90 Hz source
after save it as an avi animation
"""
import sys
sys.path.append('../source/python');

import numpy as np
import PIL
from WaveFd.Imp2DLuWave import Imp2DLuSparseWave
from WaveFd import WaveAnim
from WaveFd.Wavelet import Triangle
# loading velocity from wedge image file
filename='Wedge160x400.png'
img = PIL.Image.open(filename)
img.load()
img = img.convert('I') # gray scale, convert format
velocity  = np.asarray(img, dtype=np.float32)
velocity[:][velocity[:] == 255.0] = 3500.0
velocity[:][velocity[:] == 0.0] = 2500.0
velocity[:][velocity[:] == 9.0] = 4500.0
print np.shape(velocity)
samplerate=0.001
triangle = Triangle(100.0, samplerate)
wd2d = Imp2DLuSparseWave(400, 160, 10, samplerate, velocity, 225, 0, 1000, wavelet=triangle)
movie= wd2d.Simulate()
del wd2d # just cleaning up
movie = np.save('Imp2DWedge_400x160x1000.npy', movie) # just for preservation save it also
movie = np.load('Imp2DWedge_400x160x1000.npy')
WaveAnim.Wave2DAnim(movie, 10, 0.001, velocity, vmin=-0.002, vmax=0.002, filename='Imp2DWedge_400x160x1000')