r"""
Reads a 2D-seismic-line UKOOA P190 file (first arg file name)
and print (stdout) a txt file with lines as:

[line-name] [shot-point] [X] [Y] [Zwater]

for all geometry records on input file

execution example:

python geometry_from_ukooa_p190.py 0022_GREATERBRASILSPAN_subset3b_PSDM-vel.p190 > subset3b.geometry
"""
import sys

gfile = open(sys.argv[1])

se = ''  # line
sp = ''  # shot point
for line in gfile.readlines():

    if line[:1] == 'C':
        se = line[1:12].split()[0] # remove white spaces and new lines at once
        sp = line[18:25].split()[0]
        x = line[46:55].split()[0]
        y = line[55:64].split()[0]
        wdepth = line[64:70].split()[0]
        print '%14s %7s %7s %7s %7s' % (se, sp, x, y, wdepth)

