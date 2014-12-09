r"""
Reads a 2D-seismic-line VELF file (first arg file name)
and print (stdout) a txt file with lines as:

[line-name] [shot-point] [depth/time] [velocity]

for all velocity records on input file

execution example:

python linespzv_from_velf.py 0022_GREATERBRASILSPAN_subset3b_PSDM-vel.txt > subset3b_PSDM-spzv
"""
import sys

vfile = open(sys.argv[1])

se = ''  # line
sp = ''  # shot point
for line in vfile.readlines():

    if line[:4] == 'VELF':
        line = line[20:71]
        vels = [line[i:i+5] for i in xrange(0, len(line), 5)]
        for z, v in zip(vels[::2], vels[1::2])  :
            print '%14s %9s %9s %7s' % (se, sp, z, v)
    elif line[:4] == 'LINE':
        se = line[4:].split()[0]  # remove white spaces and new lines at once
    elif line[:4] == 'SPNT':
        sp = line[4:].split()[0]

