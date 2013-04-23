import sys
import numpy
import pylab
from pylab import cm 
sys.path.append('../../python');
from WaveFd import Exp2DWave
from WaveFd import Triangle

velocity = numpy.zeros([100, 100]) + 2500
wavelet = Triangle(fc=40, dt=0.0005)
w2w = Exp2DWave.Exp2DWave(100, 100, 5, 0.0005, velocity, 50, 50, 1000, 1, wavelet)
movie = w2w.Simulate()
numpy.shape(movie)

pylab.imshow(movie[808], cmap=cm.Greys_r, vmin=-0.005, vmax=0.005)
