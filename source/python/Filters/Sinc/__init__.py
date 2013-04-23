"""
Independent module. Compact version of Filters.

To play use SincLowPass()
to create the filter kernel 
and ConvFft() to apply it

Eg.

lowpass.ConvFft(Sd, lowpass.SincLowPass(1001, 0.5, 0.1) )

filters the Sd array using a LowPass filter of 1001 points
with a cut-off frequency of 0.5 Hz and sample rate of 0.1 s

You can also use a hanning taper window 
before applying the filter

Eg.
fkernel = lowpass.SincLowPass(1001, 0.5, 0.1) 
fkernel_tapered = Fkernel*WindowHann(np.size(fkernel))
...

"""

from LowPass import *