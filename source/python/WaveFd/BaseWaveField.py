r"""

*Base Classes*

"""

import sys, copy
import numpy as np
from Wavelet import RickerSource
import time

class BaseWave2DField(object):
    r"""
    Base class for 2d wave field, considers:

    * constant density : provides just a velocity field
    * dx=dz=ds discretization in x and z equal

    Provides:

    * a base wavelet source : Wavelet
    * a 2d velocity field   : Vel
    * current time field    : Ucurrent
    * future  time field    : Ufuture        

    Must implement:

    * SolveNextTime(self).
    """        
    def __init__(self,
                 nx,
                 nz,
                 ds,
                 dt,
                 velocity,
                 sx,
                 sz,
                 maxiter,
                 wavelet,
                 nrec=1
                 ):
        r"""
        Initialize a new wave equation field

        * nx       : number of discretization in x
        * nz       : number of discretization in z
        * ds       : dx=dz=ds grid spacing 
        * dt       : time step - e.g. seconds
        * velocity : 2d velocity distribution
        * sx/sz    : source wavelet position in indexes (i, k)
        * maxiter  : total iterations        
        * wavelet  : source wavelet function applied at the position (sx, sz)
        * nrec     : recording interval 1 equals time step 
        """
        self.Ds = ds
        self.Dt = dt
        self.Nrec = nrec
        self.MaxIter = maxiter
        self.Si = sx
        self.Sk = sz
        self.Wavelet = wavelet        
        self.Nx= nx
        self.Nz= nz
        
        # amplitude or strain values, for each (x, z) point 
        self.Ucurrent = np.zeros([self.Nz, self.Nx]) # current time
        self.Ufuture = np.zeros([self.Nz, self.Nx]) # future time
        # nx is like collums (i)
        # nz is like lines (k)
        # nt is like 4th dimension (t)
        # eg for firt time grid, 2nd line and 3rd colum
        # grid[0][1][2] = 0.0
        # grid[t][k][i] = 0.0        
        
        # velocity field for each (x, z) point
        # velocity doesnt vary with time
        self.Vel = np.zeros([self.Nz, self.Nx])
        # if a matrix of velocity is passed fills it        
        # setts up the velocity field
        if(type(velocity) is np.ndarray):
            if(np.shape(velocity) == (self.Nz, self.Nx) ):
                self.Vel = velocity
        # if not put a constant velocity
        else:
            self.Vel[:][:] = velocity

        # get global max and min        
        self.vmin = self.Vel[0][0]
        self.vmax = self.Vel[0][0]
        for array in self.Vel:
            vmin = array.min()
            vmax = array.max()
            if(vmin < self.vmin):
                self.vmin = vmin  
            if(vmax > self.vmax):
                self.vmax = vmax

        if not isinstance(wavelet, RickerSource):
            raise Exception(wavelet, "not a supported wavelet class")          
            # using the principle of planar waves 
            # to set the source wavelet frequency               
            # minvelocity = self.Vel[0][0]
            # for array in self.Vel:
            #     vmin = array.min()
            #     if(minvelocity > vmin):
            #         minvelocity = vmin                
            # #default wavelet triangular
            # self.Fw = minvelocity/(2*self.Ds);    
            # self.Wavelet=10*Triangle(self.Fw,self.Dt)          
            # frequency of source is not that simple to set in a base class 
        
        #print "total wavelet time : %.1f miliseconds" % (self.Dt*np.size(self.Wavelet)*1000)

    # def Source(self, tstep):
    #     r"""
    #     Set the source perturbation value at is position (Si, Sk) at a given 
    #     iteration time step.

    #     * tstep   :  given time step        
    #     """
    #     # # no perturbation after source finished
    #     # if( tstep > np.size(self.Wavelet)-1):
    #     #     return
                    
    #     # set at the current time field   
    #     self.Ufuture[self.Sk][self.Si] += self.Wavelet(tstep*self.Dt)
        
        
    def SolveNextTime(self):
        r"""
        Should implement the update of Ufuture, Ucurrent and
        any other necessary logic of your ag.
        Should create  

        * tstep variable : time step indexing, supporting the
          following logic, inexistence of variable will restart
          the process.

         ''' 
            try:
                self.tstep += 1
            except :
                self.tstep = 0
         '''       
        """
        pass
    
#    def Rewind(self):        
#        """
#        Put the simulation wave field to the begging stage        
#        """
#        try :
#            del self.tstep        
#        except:
#            pass
#        
#        self.Ucurrent[:][:] = 0.0 
#        self.Ufuture[:][:] = 0.0
        
                    
    def Simulate(self):
        r"""
        Loop through time steps solving the field.
        Uses method:. SolveNextTime().
        Prints a progress status and time spent at the end.

        * returns: A 'movie' matrix dimensions [MaxIter][Nz][Nx]            
        """
        snapiter=0
        initial = time.clock()                
        movie = np.zeros([int(self.MaxIter/self.Nrec), self.Nz, self.Nx])
        # for little problems with the wavelet put initialize as 1?
        
        for i in range(1, self.MaxIter, 1):
            self.SolveNextTime()
            if(i%self.Nrec==0): # every n'th Nrec
                movie[snapiter][:][:] = self.Ucurrent[:][:]
                snapiter+=1
            sys.stdout.write("\r progressing .. %.1f%%" %(100.0*float(i)/self.MaxIter))
            sys.stdout.flush()        
        sys.stdout.write(" done! \n")
        final = time.clock()
        sys.stdout.write("solving time (s) %.1f" %(final-initial))
                
        return movie        
       
       
    def EstimateTime(self):
        r"""
        Estimate time based on 100 interactions
        """       
        # a shallow copy will work
        clone = copy.copy(self) 
        initial = time.clock()  
        count=0          
        while (count < 100):
            clone.SolveNextTime()
            count += 1
        final = time.clock()
        timeperstep = (final-initial)/100.0        
        print "Time per step (s)", timeperstep
        print "Estimated time (s)", self.MaxIter*timeperstep        
        # let it free for the gc collector
        del clone    
        
    


class BaseWave1DField(object):
    r"""
    Base class for 1d wave field, considers:

    * constant density : provides just a velocity field
    * dx=ds discretization in x

    Provides:

    * a base wavelet source : Wavelet
    * a 1d velocity field   : Vel
    * current time field    : Ucurrent
    * future  time field    : Ufuture        

    Must implement:

    * SolveNextTime(self).
    """        
    def __init__(self,
                 nx,
                 ds,
                 dt,
                 velocity,
                 sx,
                 maxiter,
                 nrec=1,
                 wavelet=None,
                 ):
        r"""
        Initialize a new wave equation field

        * nx       : number of discretization in x
        * ds       : dx=dz=ds grid spacing 
        * dt       : time step - e.g. seconds
        * velocity : 1d velocity distribution
        * sx       : source wavelet position in indexes i
        * maxiter  : total iterations
        * nrec     : recording interval 1 equals time step 
        * wavelet  : source wavelet function applied at the position sx
          must have sample rate equal to dt
        """
        self.Ds = ds
        self.Dt = dt
        self.Nrec = nrec
        self.MaxIter = maxiter
        self.Si = sx
        self.Wavelet = wavelet        
        self.Nx= nx
        
        # amplitude or strain values, for each (x) point 
        self.Ucurrent = np.zeros(self.Nx) # current time
        self.Ufuture = np.zeros(self.Nx) # future time
        # velocity field for each x point
        # velocity doesnt vary with time
        self.Vel = np.zeros(self.Nx)
        # if a array of velocity is passed fills it        
        # setts up the velocity field
        if(type(velocity) is np.ndarray):
            if(np.shape(velocity) == (self.Nx) ):
                self.Vel = velocity
        # if not put a constant velocity
        else:
            self.Vel[:][:] = velocity


        if(wavelet == None):         
            # using the principle of planar waves 
            # to set the source wavelet frequency               
            minvelocity = self.Vel[0][0]
            for array in self.Vel:
                vmin = array.min()
                if(minvelocity > vmin):
                    minvelocity = vmin                
            #default wavelet triangular
            self.Fw = minvelocity/(2*self.Ds);    
            self.Wavelet=10*Triangle(self.Fw, self.Dt)           
        
        print "total wavelet time : %.1f miliseconds" % (self.Dt*np.size(self.Wavelet)*1000)

    def Source(self, tstep):
        r"""
        Set the source perturbation value at is position (Si) at a given 
        iteration time step.

        * tstep   :  given time step        
        """
        # no perturbation after source finished
        if( tstep > np.size(self.Wavelet)-1):
            return
                    
        # set at the current time field   
        self.Ucurrent[self.Si] = self.Wavelet[tstep]
        
        
    def SolveNextTime(self):
        r"""
        Should implement the update of Ufuture, Ucurrent and
        any other necessary logic of your ag.

        Should create:
 
        * tstep variable : time step indexing, supporting the
          following logic, inexistence of variable will restart
          the process.
         ''' 
            try:
                self.tstep += 1
            except :
                self.tstep = 0
         '''       
        """
        pass


    def Simulate(self):
        r"""
        Loop through time steps solving the field.
        Uses method:. SolveNextTime().
        Prints a progress status and time spent at the end.

        * returns: A 'movie' matrix dimensions [MaxIter][Nx]            
        """
        snapiter=0
        initial = time.clock()                
        movie = np.zeros([int(self.MaxIter/self.Nrec), self.Nx])


        for i in range(1, self.MaxIter, 1):
            self.SolveNextTime()
            if(i%self.Nrec==0): # every n'th Nrec
                movie[snapiter] = self.Ucurrent
                snapiter+=1
            sys.stdout.write("\r progressing .. %.1f%%" %(100.0*float(i)/self.MaxIter))
            sys.stdout.flush()        
        sys.stdout.write(" done! \n")
        final = time.clock()
        sys.stdout.write("solving time (s) %.1f" %(final-initial))
                
        return movie        
       
       
    def EstimateTime(self):
        r"""
        Estimate time based on 50 interactions
        """       
        # a shallow copy will work
        clone = copy.copy(self) 
        initial = time.clock()    
        count=0        
        while (count < 50):
            clone.SolveNextTime()
            count +=1
        final = time.clock()
        timeperstep = (final-initial)/50.0        
        print "Time per step (s)", timeperstep
        print "Estimated time (s)", self.MaxIter*timeperstep        
        # let it free for the gc collector
        del clone    
        
    
