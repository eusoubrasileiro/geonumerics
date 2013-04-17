#!/usr/bin
# import filters.py

#backend = 'gtk'
#import matplotlib
#matplotlib.use(backend)

import pylab as py
import numpy as np
import time

"""
Explicit wave equation (acoustic) , finite differences
3rd order centered in space
2 order backward in time
convergence limitations
"""

class W2DWExp:
    def __init__(self, Nx=100, Nz=100, Dx=0.5, Dz=0.5, Dt=0.004):
        """
        initialize a new wave equation field,
        for solving with finite diferences method
        Nx number of discretization in x  - ground dimension (e.g. meters)
        Nz number of discretization in z
        Dx grid spacing in x
        Dz grid spacing in z
        Dt time step (e.g. seconds)
        """
        self.Nx=Nx
        self.Nz=Nz
        self.Dx=Dx
        self.Dz=Dz
        self.Dt=Dt
        # 3rd order on space, centered
        # due 3rd order finite diferences ... N+1 + order
        # real dimensions are N+1+3 = N+4
        self.nx = Nx+4
        self.nz = Nz+4
        # 2nd order time, backward
        # so we need plus order+1 grids
        
        # amplitude or strain values, for each (x, z) point
        # must follow the order nt, nz, nx
        self.U0 = np.zeros([self.nz, self.nx]) # t-1
        self.U1 = np.zeros([self.nz, self.nx]) # t
        self.U2 = np.zeros([self.nz, self.nx]) # t+1

        # nx is like collums (i)
        # nz is like lines (k)
        # eg for firt time grid, 2nd line and 3rd colum
        # self.U0[1][2] = 0.0
        # self.U0[k][i] = 0.0
        # velocity field for each (x, z) point
        # velocity doesnt vary with time
        self.Gvel = np.zeros([self.nz, self.nx])

    def SetVel(self,Velocity=2000):
        """
        sets the velocy field
        """
        # if a matrix of velocity is passed fills it
        if(type(Velocity) is np.ndarray):
            if(np.shape(Velocity) == (self.nz, self.nx) ):
                self.Gvel = Velocity
        # if not put a constant velocity
        else:
             self.Gvel[:][:] = Velocity

    def NextTimeIK(self, i, k):
        """
        Calculates the Laplacian L at the (i, k) position
        considering u(i, k) the displacement fielt at time t
        L(i,k) = d2u + d2u
                 dx2   dz2
        centered differences 3rd order
        """

        # __Avoiding problematic positions!!!
        # samples needed around the (i, k) position
        # depends on the order of finite diferences 
        # finite diference in space
        # 3rd order = 2 samples before and after
        if(i+2 > self.nx or i-2 < 0):
            return
        if(k+2 > self.nz or k-2 < 0):
            return
            
        # in time t
        u = self.U1

        dx = self.Dx
        d2u_dx2 = -(u[k][i+2]-2*u[k][i]+u[k][i-2])/12.0
        d2u_dx2 += 4.0*(u[k][i+1]-2*u[k][i]+u[k][i-1])/3.0
        d2u_dx2 /= dx

        dz = self.Dz
        d2u_dz2 = -(u[k+2][i]-2*u[k][i]+u[k-2][i])/12.0
        d2u_dz2 += 4.0*(u[k+1][i]-2*u[k][i]+u[k-1][i])/3.0
        d2u_dz2 /= dz

        LaplacianIK = d2u_dx2 + d2u_dz2;

        if(i == self.Si and k == self.Sk and self.t < np.size(self.SourceWavelet)):
            LaplacianIK -= self.SourceWavelet[self.t]

        U1IK = LaplacianIK*(self.Gvel[k][i]*self.Dt)**2
        U1IK += 2*self.U1[k][i]-self.U0[k][i]

        return U1IK

    def SetSource(self, i, k, Wavelet):
        self.Si = i
        self.Sk = k
        self.SourceWavelet = Wavelet        
        self.t = 1

        return


    def NextTime(self):
        """
        same as before but...
        for the entire displacement field,
        after the process U1 will be
        replaced by U2 ... and so on
        """

        # next time U2
        U = self.U2

        # ignore external part of the grid, zero fixed
        # 3rd order fixed
        for i in range(2, self.Nx-2, 1):
            for k in range(2, self.Nz-2, 1):
                U[k][i]=self.NextTimeIK(i,k)

        # make the update in the time stack
        self.U0 = self.U1
        self.U1 = U        


    def MoveForward(self):
        """
        t must be set as 1, first time step
        """
        
        self.NextTime()        
    
        self.t = self.t + 1

def main():
    field = W2DWExp(10,10,Dt=0.1)
    # 10 m/s
    field.SetVel(0.04)
    # initial condition at t
    # t is 2
    field.UTime[0][5][5]=0.0
    field.UTime[1][5][5]=1.0
    py.ion()
    img = py.imshow(field.UTime[1])
    py.show()

    for i in range(50):
        field.next_time_all()
        img.set_data(field.UTime[1])
        py.draw()
#        time.sleep(0.2)
    return field.UTime[1]



if __name__ == '__main__':
    main()
