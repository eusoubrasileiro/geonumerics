#!/usr/bin/python
import numpy
from scipy.interpolate import InterpolatedUnivariateSpline
from obspy.segy.core import readSEGY
import os

def SEGYcoordinatesfix(segyfilename):
    section = readSEGY(segyfilename, unpack_trace_headers=True)
    # getting just not repeated coordinate values sx, sy
    ntr = len(section) # number of traces
    sx = numpy.zeros(1)
    sy = numpy.zeros(1)
    trc = numpy.zeros(1) # trace index of not duplicated traces
    cdpx = numpy.zeros(ntr) # original x coordinate
    cdpy = numpy.zeros(ntr) # original y coordinate
    # bytes (181, 185) (cdpx, cdpy) (first point allways in)
    cdpx[0] = section[0].stats.segy.trace_header.x_coordinate_of_ensemble_position_of_this_trace
    cdpy[0] = section[0].stats.segy.trace_header.y_coordinate_of_ensemble_position_of_this_trace
    sx[0] = cdpx[0]
    sy[0] = cdpy[0]
    trc[0] = 0
    for i in numpy.arange(1, ntr): # get just the not duplicated coordinates
        cdpx[i] = section[i].stats.segy.trace_header.x_coordinate_of_ensemble_position_of_this_trace
        cdpy[i] = section[i].stats.segy.trace_header.y_coordinate_of_ensemble_position_of_this_trace
        if (cdpx[i] != cdpx[i-1]) or (cdpy[i] != cdpy[i-1]):  # just when (x, y) == (x, y) ignore
            sx = numpy.append(sx, cdpx[i])
            sy = numpy.append(sy, cdpy[i])
            trc = numpy.append(trc, i)
    #trc (not duplicated indexes = x)
    #sx, sy not duplicated coordinates
    flinearsx = InterpolatedUnivariateSpline(trc, sx, bbox=[-3, ntr+2], k=1) # linear iterp function on xcoordinate ; x is trace index
    flinearsy = InterpolatedUnivariateSpline(trc, sy, bbox=[-3, ntr+2], k=1) # linear iterp function on ycoordinate ; x is trace index
    # (to enable linear extrapolation that interp1 doesn't do) spline=linear iterp function case where spline degree k=1
    # uses limits of extrapolation +3 traces before and after
    for trace_index in numpy.arange(0, ntr, 1): # interpolate for all trace indexes, changing the trace headers on bytes (73, 77)
        section[trace_index].stats.segy.trace_header.source_coordinate_x = int(flinearsx(trace_index))
        section[trace_index].stats.segy.trace_header.source_coordinate_y = int(flinearsy(trace_index))
    fileName, fileExtension = os.path.splitext(segyfilename)
    section.write(fileName+'fixed.segy', format='SEGY')

# fix the entire list of files terminated with *.sgy in the directory where the script runs
import glob

print "This script fixes all coordinates duplicated in the trace header."
print "It read x, y from bytes (181, 185) and writes on bytes (73, 77)"
print "It will run for of all segy's in the current folder"
print ": extension is *.sgy"

files = glob.glob('./*.sgy')
i = 1
for asegy in files:
    print "file ", asegy, " ",  i, " of ", len(files)
    i += 1
    SEGYcoordinatesfix(asegy)

