import sys

import os

import salome

import numpy

from salome.geom import geomtools


def write_points( geompy, points ):

    directory = "lattice"
    
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    file = open(directory + "/points", 'w')


    # Number of points
    
    file.write( str(len(points)) )
    file.write( "\n" )


    for point in points:

        coord = geompy.PointCoordinates(point)

        file.write( "%.0f " % coord[0] )
        file.write( "%.0f " % coord[1] )
        file.write( "%.0f\n" % coord[2] )


    file.close()
