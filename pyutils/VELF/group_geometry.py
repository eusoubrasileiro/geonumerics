__author__ = 'andre'
r"""
group many sp geometry files in an unique file
printing a new collum with the name of the file
that should be the line name
"""
import glob

files = glob.glob('*.spgeo')
i = 1
for asegy in files:
    fs = open(asegy, 'r')
    for line in fs.readlines():
        print '%9s     %s' % (asegy.split('.spgeo')[0], line.strip())
    fs.close()
    i += 1