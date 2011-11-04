#!/usr/bin
# import filters.py

#backend = 'gtk'
#import matplotlib
#matplotlib.use(backend)

import pylab as py
import numpy as np
import time

__doc__ = """


"""

class wave_equation:
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
        self.fospace = 3
        self.nx = Nx+(self.fospace+1)
        self.nz = Nz+(self.fospace+1)
        # 2nd order time, backward
        # so we need plus order+1 grids
        self.fotime=2
        self.nt = self.fotime+1
        
        # amplitude or strain values, for each (x, z) point
        # must follow the order nt, nz, nx
        self.UTime = np.array(np.zeros([self.nt, self.nz, self.nx]), dtype=float)
        # nx is like collums (i)
        # nz is like lines (k)
        # nt is like 4th dimension (t)
        # eg for firt time grid, 2nd line and 3rd colum
        # grid[0][1][2] = 0.0
        # grid[t][k][i] = 0.0
        # velocity field for each (x, z) point
        # velocity doesnt vary with time
        self.Gvel = np.array(np.zeros([self.nz, self.nx]), dtype=float)

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


    def laplacian(self, i, k):
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
        maxn = self.fospace-1
        if(i+maxn > self.nx and i-maxn<0):
            return
        if(k+maxn > self.nt and k-maxn<0):
            return
            
        # in time
        # As t is in [0, 1, 2] (2nd order)
        # time t in this case is
        # UTime[self.fotime-1] = UTime[2-1] = UTime[1]
        u = self.UTime[1]

        dx = self.Dx
        d2u_dx2 = -(u[k][i+2]-2*u[k][i]+u[k][i-2])/12.0
        d2u_dx2 += 4.0*(u[k][i+1]-2*u[k][i]+u[k][i-1])/3.0
        d2u_dx2 /= dx

        dz = self.Dz
        d2u_dz2 = -(u[k+2][i]-2*u[k][i]+u[k-2][i])/12.0
        d2u_dz2 += 4.0*(u[k+1][i]-2*u[k][i]+u[k-1][i])/3.0
        d2u_dz2 /= dz

        return d2u_dx2 + d2u_dz2;


    def next_time(self, i, k):
        """
        Calculate the time t from:
        a) the previous (t-2, t-1, t)
        b) laplacian at time t

        """
        # __Avoiding problematic positions!!!
        # samples needed around the (i, k) position
        # depends on the order of finite diferences 
        # finite diference in space
        # 3rd order = 2 samples before and after
        maxn = self.fospace-1
        if(i+maxn > self.nx and i-maxn<0):
            return
        if(k+maxn > self.nt and k-maxn<0):
            return

        # As t is in [  0,   1,  2] (2nd order)
        #            [t-2, t-1,  t]
        u = self.UTime
        dt = self.Dt

        ut = self.laplacian(i,k)*(self.Gvel[k][i]*dt)**2
        ut += -2*u[1][k][i]+u[0][k][i]

        return ut



    def next_time_all(self):
        """
        same as before but...
        for the entire displacement field,
        after the process Utime[t+1] will be
        replaced by Utime[t+2] ... and so on
        """


        # in time
        # As t is in [0, 1, 2] (2nd order)
        # time t in this case is
        # UTime[self.fotime-1] = UTime[2-1] = UTime[1]
        u = self.UTime[2]

        # ignore external part of the grid
        # 3rd order, self.order-1
        for i in range(self.fospace-1, self.fospace-1+self.Nx, 1):
            for k in range(self.fospace-1, self.fospace-1+self.Nz,1):
                u[k][i]=self.next_time(i,k)

        # make the update in the time stack
        # before [t-2, t-1,  t]
        # after  [t-1, t,  t+1]
        # so t-2 receive t-1 and etc...        

        u = self.UTime # points to general pointer again
        u[0] = u[1]
        u[1] = u[2]
        
        return

    def initial_condition():
        pass


def main():
    field = wave_equation(10,10,Dt=0.1)
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
r#        time.sleep(0.2)
    return field.UTime[1]



if __name__ == '__main__':
    main()
