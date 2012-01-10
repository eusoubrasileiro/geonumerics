'''
Created on Dec 29, 2011

@author: andre
'''
import unittest
from Filters import samplesignals


class Test(unittest.TestCase):

    def setUo(self):
        self.signal = samplesignals.Periodic(200, 0.01)
        
    def testName(self):
        pass


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
