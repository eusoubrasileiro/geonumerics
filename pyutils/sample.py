#!/usr/bin/python
import sys
import numpy

for i in numpy.arange(1,len(sys.argv)):
	print sys.argv[i].split('.')[0]+'B.segy'

