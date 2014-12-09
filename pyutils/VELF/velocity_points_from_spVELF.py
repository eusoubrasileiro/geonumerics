r"""
From a line-sp velocity file `linespzv_from_velf.py` (first arg file name)
and a geometry-file from 2D-seismic-line obtained with Seisee
from the trace header information (sp-x-y) and added a line name column (second arg file name)

print (stdout) a txt file with lines as:

[line-name] [shot-point] [X] [Y] [depth/time] [velocity]

for velocity records on input file

execution example:

python velocity_points_from_spVELF.py subset3bn5-PSDM-spzv jacuipe-gxt.geox > jacuipe_vels.points

"""
import sys
from scipy.interpolate import InterpolatedUnivariateSpline
import pandas

vfile = sys.argv[1]
gfile = sys.argv[2]

geometry_table = pandas.read_table(gfile, header=None, names=['line', 'sp', 'x', 'y'], skipinitialspace=True, sep=' ')
linespvel_table = pandas.read_table(vfile, header=None, names=['line', 'sp', 'z', 'vel'], skipinitialspace=True, sep=' ')

# build interpolators for every line on geometry file
groupedbylines = geometry_table.groupby('line')
interpolators = []
for linename, group in groupedbylines:
    x = group['x']
    y = group['y']
    sp = group['sp']
    # linear iterp function on (x,y) coordinates as function of sp (extrapolation at maximum)
    flinearx = InterpolatedUnivariateSpline(sp, x, bbox=[min(sp)-min(sp), 2*max(sp)], k=1)
    flineary = InterpolatedUnivariateSpline(sp, y, bbox=[min(sp)-min(sp), 2*max(sp)], k=1)
    interpolators.append({'linename': linename, 'fx': flinearx, 'fy': flineary})

# interpolate the velocity file based on the sp coordinates
for i in xrange(len(linespvel_table)):
    linename = linespvel_table['line'][i]
    sp = linespvel_table['sp'][i]
    z = linespvel_table['z'][i]
    vel = linespvel_table['vel'][i]

    try:
        funcs = filter(lambda x: x['linename'] == linename, interpolators)[0]
    except:
        continue
        sys.stderr.write("Error "+str(linename)+" "+str(sp)+" "+str(z)+" "+str(vel)+"\n")
        # when comes to here might be lines without interpolation functions

    x = funcs['fx'](sp)
    y = funcs['fy'](sp)

    print '%14s %9d %11.1f %11.1f %9d %9d'%(linename, sp, x, y, z, vel)