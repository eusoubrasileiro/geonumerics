# xcorr function original from mathlab is equal a normal correlation
# padding the small vector with 00's
# the xcorr function from pylab do not allow vector being of diferrent size

import numpy, sys, os
from pylab import arange
import pylab
import pyplot

x = [0, 0, 1, 2, 3, 3, 3, 3, 2, 1, 0, 0, 0, 0, 0, 0, 0];
x_2 = [0, 0] + x; # x time shift 2 samples

pylab.plot(x);
pylab.plot(x_2);
pylab.ylim(ymax=5);
pylab.figure();
x = x + [0, 0] # add two samples on x, so now they have the same size

pylab.xcorr(x, x_2);
# so  x_2 tem que andar -2 amostras pra ficar em maxima correlacao com x
# fazendo o mesmo with correlate que e uma convolucao com o segundo invertido

corr = pylab.correlate(x, x_2, mode="full");
lags = pylab.arange(-len(x)+1, len(x), 1); # the lags, cobrem 2N-1, entao len-1 ate len-1, incluindo o 0
pylab.figure();
pylab.plot(lags, corr);