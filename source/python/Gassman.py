import numpy, sys, os
import pylab
from numpy.random import rand
import time


def RealKsaturated(Vp, Vs, constBulkDensity):
    """
    Calculate real K modules from Vp, Vs and bulk density of the rock
    Vp and Vs must be array's. (Remember units!!)
    """
    
    return constBulkDensity*(Vp**2-4*Vs**2/3)


# to make fluid substituition, first we have to define:
# 1) porosity of rock
# 2) properties of the fluids (K_fl, Density_fl)
# 3) bulk modulus of mineral matrix K0
# 4) bulk modulus of porous rock frame K*
