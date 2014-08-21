"""
2D line split spread seismic acquisition cmp pattern chart

* shots are over station positions
not stack array

* station number is increasing

"""

import numpy as np
import pylab as py

dstation = 1  # station interval (change just if wanting to use any other unit e. g. meters)
# all units (station units)
dsource = 2  # source displacement
halfspread = 20  # half spread number of channels }
linenumber = 0  # line number prefix for stations
nstations = 50  # number of stations on line
laststation = nstations*dstation + linenumber
stations = np.linspace(linenumber, laststation, nstations)  # station numbers in line
nshots = int(nstations/dsource)  # number of shots
cmps = np.ma.zeros((nshots, nstations))  # cmp coverage matrix (shots, cmp hits per shot)
for ep in stations[::dsource]:
    begin = ep - halfspread
    end = ep + halfspread
    if begin < linenumber:
        begin = linenumber
    if end > laststation:
        end = laststation
    epcmps = (np.ma.masked_outside(stations, begin, end) - ep)/2. + ep  # mid points
    cmps[(ep-linenumber)/dsource] = epcmps[:]

# make chart
py.figure()
for i in range(nshots):
    cmp = cmps[i].compressed()
    py.scatter(i*dsource, i*dsource+1, c='k', marker='*')
    py.scatter(cmp, np.zeros(cmp.size)+i*dsource, c='b', marker='+')
    py.xlabel("stations")
    py.ylabel("shot number")

py.show()