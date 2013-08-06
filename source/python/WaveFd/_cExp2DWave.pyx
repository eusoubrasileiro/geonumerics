import numpy
import sys
from libc.math cimport exp
# Import Cython definitions for numpy
cimport numpy
cimport cython

DTYPE = numpy.double
ctypedef numpy.double_t DTYPE_T

# it's fundamental to use ctypes (cdef) indexes otherwise no performance boost in loops
# remove if's inside loop to better performance

@cython.boundscheck(False)
@cython.wraparound(False)
def SolveUfuture(
    numpy.ndarray[DTYPE_T, ndim=2] Uf not None, 
    numpy.ndarray[DTYPE_T, ndim=2] Uc not None, 
    numpy.ndarray[DTYPE_T, ndim=2] Up not None, 
    numpy.ndarray[DTYPE_T, ndim=2] V not None, 
    int Nz, 
    int Nx, 
    double Ds, 
    double Dt):

    cdef int k, j
    cdef double u0k, u1k, u3k, u4k, uj0, uj1, uj3, uj4, ujk

    for k in range(Nz):
        for j in range(Nx):  
            u0k=u1k=u3k=u4k=0.0
            uj0=uj1=uj3=uj4=0.0
            ujk = Uc[k,j]

            if(j-2 > -1):
                u0k = Uc[k,j-2]   
            if(j-1 > -1):
                u1k = Uc[k,j-1]
            if(j+1 < Nx):
                u3k = Uc[k,j+1]            
            if(j+2 < Nx):
                u4k = Uc[k,j+2]
            if(k-2 > -1):
                uj0 = Uc[k-2,j]
            if(k-1 > -1):
                uj1 = Uc[k-1,j]
            if(k+1 < Nz):
                uj3 = Uc[k+1,j]            
            if(k+2 < Nz):
                uj4 = Uc[k+2,j]       
            
            Uf[k,j] = ((1.0/12.0)*(16.0*(u1k+u3k+uj1+uj3)-(u4k+uj0+u0k+uj4)-60.0*ujk)*(Dt*V[k,j]/Ds)**2)
            Uf[k,j] += 2*Uc[k,j]-Up[k,j]

