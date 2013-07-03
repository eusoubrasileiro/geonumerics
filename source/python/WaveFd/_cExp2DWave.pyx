import numpy

from libc.math cimport exp
# Import Cython definitions for numpy
cimport numpy
cimport cython

DTYPE = numpy.float
ctypedef numpy.float_t DTYPE_T

@cython.boundscheck(False)
@cython.wraparound(False)
def SolveUfuture(
	numpy.ndarray[DTYPE_T, ndim=2] Ufuture not None, 
	numpy.ndarray[DTYPE_T, ndim=2] Ucurrent not None, 
	numpy.ndarray[DTYPE_T, ndim=2] Uprevious not None, 
	numpy.ndarray[DTYPE_T, ndim=2] V not None, 
	unsigned int Nz, 
	unsigned int Nx, 
	double Ds, 
	double Dt):

	for k in range(Nz):
		for j in range(Nx):
			# u0k u1k*ujk*u3k u4k     
			# Boundary fixed 0 outside        
			u0k=u1k=u3k=u4k=0.0
			uj0=uj1=uj3=uj4=0.0
			ujk = Ucurrent[k,j]      

			if(j-2 > -1):
				u0k = Ucurrent[k,j-2]	
			if(j-1 > -1):
				u1k = Ucurrent[k,j-1]
			if(j+1 < Nx):
				u3k = Ucurrent[k,j+1]            
			if(j+2 < Nx):
				u4k = Ucurrent[k,j+2]
			if(k-2 > -1):
				uj0 = Ucurrent[k-2,j]
			if(k-1 > -1):
				uj1 = Ucurrent[k-1,j]
			if(k+1 < Nz):
				uj3 = Ucurrent[k+1,j]            
			if(k+2 < Nz):
				uj4 = Ucurrent[k+2,j]		
			
			Ufuture[k,j] = (1.0/12)*(-u0k+16*u1k+16*u3k-u4k-uj0+16*uj1+16*uj3-uj4-60*ujk)*(Dt*V[k,j]/Ds)**2 +2*ujk-Uprevious[k,j]