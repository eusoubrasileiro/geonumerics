r"""

Manipulation utilities and others

"""
from PIL import Image
import numpy as np


# weights for centered finite diferences second derivative 
centered2ndfiniteweight6th = np.array([1.0/90, -3.0/20, 3.0/2, -49.0/18, 3.0/2, -3.0/20, 1.0/90])
centered2ndfiniteweight4th = np.array([ -1.0/12, 4.0/3, -5.0/2, 4.0/3, -1.0/12])
centered2ndfiniteweight2nd = np.array([ 1.0, -2.0, 1.0 ]) 

def Traces2DFromSimulation(simulation, ds, s0, ntrac):
    """
    From a simulation movie matrix creates a matrix
    of traces equally spaced...

    * simulation  : simulation matrix[steps][nz][nx]
    * ds          : cell spacing x/z 
    * s0          : first trace positon in world coordinates
    * ntrac       : number of traces starting at s0
    
    Returns:
    
    * (extent, traces[t](i)) : limits for using imshow and traces

    """
    nt = np.shape(simulation)[0] 
    s0 = s0/ds
    traces = np.zeros([nt, ntrac])
    for i in range(ntrac):
        si= s0 + i 
        sk= 0    
        for t in range(nt):
            traces[t][i] = simulation[t][sk][si]
    xmin, xmax, ymin, ymax = s0*ds, s0*ds+ntrac*ds, 0, nt
    extent = xmin, xmax, ymax, ymin

    return extent, traces

def LoadPicture( filename ):    
    img = Image.open(filename)
    img = img.convert('L') # gray scale, convert format
    data  = np.asarray(img, dtype=np.float32)
    return data

