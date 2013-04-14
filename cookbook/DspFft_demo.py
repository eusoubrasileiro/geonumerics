#-------------------------------------------------------------------------------
# Name:        Dft vs Fft c code demos
# Purpose:     some examples using the Filters.py library
#
# Author:      andre
#-------------------------------------------------------------------------------
#!/usr/bin/env python

import sys
sys.path.append('..\\source\\python');
sys.path.append('..\\source\interfaces'); #where is the dsprocessing.pyd (dll)
# or to not windows guys
# sys.path.append(os.path.join('..','log','python'))
# the above to be able to load the modules bellow
import samplesignals
import filters
import pylab
import numpy
from numpy.random import rand
import time
import dsprocessing

def main():
    pylab.ion();
    ind = [0,];
    ldft = [0,];
    lfft = [0,];
    lpfft = [0,]

    # plot a graph Dft vs Fft, lists just support size until 2**9
    for i in range(1, 9, 1):
        t_before = time.clock();
        dsprocessing.dspDft(rand(2**i).tolist());
        dt = time.clock() - t_before;
        ldft.append(dt);
        print ("dft ", 2**i, dt);
        #pylab.plot([2**i,], [time.clock()-t_before,]);
        t_before = time.clock();
        dsprocessing.dspFft(rand(2**i).tolist());
        dt = time.clock() - t_before;
        print ("fft ", 2**i, dt);
        lfft.append(dt);
        #pylab.plot([2**i,], [time.clock()-t_before,]);
        ind.append(2**i);
        # python fft just to compare
        t_before = time.clock();
        pylab.fft(rand(2**i).tolist());
        dt = time.clock() - t_before;
        lpfft.append(dt);

    pylab.plot(ind, ldft);
    pylab.plot(ind, lfft);
    pylab.plot(ind, lpfft);
    pylab.show();
    return [ind, ldft, lfft, lpfft];

if __name__ == '__main__':
    main()

