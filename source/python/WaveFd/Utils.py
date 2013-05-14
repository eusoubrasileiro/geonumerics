r"""

Manipulation utilities and others

"""
import Image
import numpy as np

def Traces2DFromSimulation(simulation, ds):
    """
    From a simulation movie matrix creates a matrix
    of traces of equally spaced...

    name  :  file name begin of the image files *.png
             also the name of the output *avi file
    gif   :  if False create an *.avi file instead of a *.gif
    fps   :  frames per second

    """

def LoadPicture( filename ):
    img = Image.open( filename )
    img.load()
    data = np.asarray( img, dtype="float32" )
    return data

