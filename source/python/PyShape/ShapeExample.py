# simple example in how to create a simple shape file
# 2 lines, coordinates of Brazil Es Basin, minimu 5 points
import shapefile
w = shapefile.Writer(shapefile.POLYGON)
w.field('FIRST_FLD','C','40') # corresponds to 'First', 40 characteres
w.field('SECOND_FLD','C','40') # corresponds to 'Polygon' 40 characteres
parts= [[[-39.8,-20.58],[-39.2,-19.69],[-38.74,-20.47], [-38.70,-20.40], [-38.65,-20.35]]]
w.line(parts)
parts = [[[-38.79,-19.70],[-38.64,-19.78]]]
w.line(parts)
w.record('First','Line') # first record corresponds to the first added shape
w.record('Second','Line') # second record corresponds to the seconds added shape
# number of records must be equal to the number of shapes
w.save('TwoLines.shp')

# more about ... autoBalance garantees we get everything balanced in booth sides
#>>> w.autoBalance = 1
#You also have the option of manually calling the balance() method each time you add a shape or a record
#to ensure the other side is up to date. When balancing is used null shapes are created on the geometry
#side or a record with a value of "NULL" for each field is created on the attribute side.