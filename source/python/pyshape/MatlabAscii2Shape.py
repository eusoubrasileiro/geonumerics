# script to create a shape file from a matlab file
# using the shape file python library http://code.google.com/p/pyshp/
# turns a matlab polyline vector format in an shape file
# file downloaded from http://www.ngdc.noaa.gov/mgg/coast/
# using matlab coast format


import shapefile
import math

def CreateShape(MatlabFile, ShapeName='Shape.shp', Type=shapefile.POLYLINE):    
    """
    usage:
    CreateShape('Matlabfile.dat', 'Sample.shp', shapefile.POINT)
    
    where the MatlabFile text must have 
    the following format
    '
    -70.034841	0
    -70.034547	0.075979
    -69.999931	0.554733
    nan nan 
    -69.220782	1.000045
    -70.034547	0.075979
    nan nan 
    '
    nan values 
    must be here in the end also
    
    The available shape file types are 
    shapefile.POINT
    shapefile.POLYLINE
    shapefile.POLYGON
    shapefile.MULTIPOINT 
    
    """
    # parse the file to an array  
    
    fp = open(MatlabFile, 'r')
    file = []
    for line in fp.readlines():
        line = line.replace('\t', ' ') #guarante split character equal space ' '
        line = line.replace(';', ' ')
        file.append([float(line.split(' ')[0]), float(line.split(' ')[1])])
    
    # create the shape and fill it
    w = shapefile.Writer(Type)
    w.field('FIRST_FLD','C','40')
    # an simple field of 40 characters
    parts = []     
    
    # pointer to the correct function
    pcreateShape = []    
    if(Type == shapefile.POINT or Type == shapefile.MULTIPOINT):
        # add each point as points
        # ignore nan values
        for pair in file:
            if(math.isnan(pair[0])):                
                continue
            w.point(pair[0], pair[1])
    
    if(Type == shapefile.POLYLINE):
        # add each chunk of points as lines
        #chunks diveded by nan values
        for pair in file:
            if(math.isnan(pair[0])):
                w.line([parts])
                parts = []
                continue
            parts.append(pair)
    
    if(Type == shapefile.POLYGON):
        # add each chunk of points as polygons
        #chunks diveded by nan values
        for pair in file:
            if(math.isnan(pair[0])):
                w.poly([parts])
                parts = []
                continue
            parts.append(pair)
    
    # guarantees autobalance between records and shapes, not working
    # w.autoBalance = 1
    # now working, my own balance
    for i in range(len(w.shapes())):
        w.record('Null')
        
    print len(w.shapes())
    print len(w.records)
    w.save(ShapeName)  


