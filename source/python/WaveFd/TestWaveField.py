'''
Created on Dec 29, 2011

@author: andre
'''
import unittest
import pylab as py
import numpy as np
from WaveFd import WaveField
import time


class Test(unittest.TestCase):

    def testName(self):
        pass

    def testExampleLayers(self):
        fd = WaveField(100,50,0.5,0.005,Fw=10)
        fd.MaxIter=20
        # 100*0.5 = 50meters
        # 20*0.5 = 10meters
        # 10 m/s
        fd.SetVel(10)
        fd.Vel[25:49][:]=3
        # initial condition at t
        # t is 2
        #field.Utime[1][1][1]=100.0
        #field.Utime[0][1][1]=0
        return fd.Loop()

def Plot():
    py.ion()
    ti = np.load('IfE'+str(1)+'.npy')
    img = py.imshow(np.reshape(ti, (50,100))) #@UnusedVariable
    
    for i in range(2,10):
        ti = np.load('IfE'+str(i)+'.npy')
        py.imshow(np.reshape(ti, (50,100)))
        py.draw()
        time.sleep(0.1)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
