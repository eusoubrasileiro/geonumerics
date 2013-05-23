r"""

4th order explicit space test question

Convergence
 
.. math::
 
     \Delta t \leq \frac{2 \Delta s}{ V \sqrt{\sum_{a=-N}^{N} (|w_a^1| + |w_a^2|)}}

Where w_a are the centered differences weights

Under limits.

"""

import numpy as np
import sys

def Triangle(fc, dt, n=None):
    r"""
    Triangle Wave one Period.
    Defined by frequency and sample rate or by size

    * fc        : maximum desired frequency
    * dt        : sample rate
    * n         : half length of triangle    
    """
    if(n==None):
        n=int(1/float(fc*dt))

    t = np.arange(0+1.0/n, 1, 1.0/n)
    y = 1-t
    y = np.append(y, 0.0)
    y_ = 1-t[::-1]
    y_ = np.insert(y_, 0, 0.0)
    
    return np.append(y_, np.append(1, y))

Nz = Nx = 20
Dt = 0.001
Ds = 10
numberiter = 200

Source = Triangle(90, 0.001)
Uprevious = np.zeros([Nz, Nx]) # previous time
Ucurrent = np.zeros([Nz, Nx]) # current time
Ufuture = np.zeros([Nz, Nx]) # future time
V = np.zeros([Nz, Nx]) 
V[:][:] = 2000.0

# additional not needed 
Simulation  = np.zeros([numberiter, Nz, Nx])

# source activation, center of grid
Uprevious[10][10] = Source[0]

for i in range(1, numberiter+1):

	# tringular source position center grid
	if(i < np.size(Source)):
		Ucurrent[10][10] = Source[i]

	for k in range(Nz):
		for j in range(Nx):
			# u0k u1k*ujk*u3k u4k     
			# Boundary fixed 0 outside        
			u0k=u1k=u3k=u4k=0.0
			uj0=uj1=uj3=uj4=0.0
			ujk = Ucurrent[k][j]      

			if(j-2 > -1):
				u0k = Ucurrent[k][j-2]	
			if(j-1 > -1):
				u1k = Ucurrent[k][j-1]
			if(j+1 < Nx):
				u3k = Ucurrent[k][j+1]            
			if(j+2 < Nx):
				u4k = Ucurrent[k][j+2]
			if(k-2 > -1):
				uj0 = Ucurrent[k-2][j]
			if(k-1 > -1):
				uj1 = Ucurrent[k-1][j]
			if(k+1 < Nz):
				uj3 = Ucurrent[k+1][j]            
			if(k+2 < Nz):
				uj4 = Ucurrent[k+2][j]

			d2u_dx2 = (-u0k+16*u1k-30*ujk+16*u3k-u4k)/12.0
			d2u_dz2 = (-uj0+16*uj1-30*ujk+16*uj3-uj4)/12.0
			Ufuture[k][j] = (d2u_dx2+d2u_dz2)*(Dt*V[k][j]/Ds)**2
			Ufuture[k][j] += 2*Ucurrent[k][j]-Uprevious[k][j]

	# make the update in the time stack
	Uprevious = Ucurrent
	Ucurrent = Ufuture
	Simulation[i-1] = Ucurrent

	sys.stdout.write("\r %d" %(i) )