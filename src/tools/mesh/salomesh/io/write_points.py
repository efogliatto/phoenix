import sys

import os

import salome

import numpy

# from salome.geom import geomtools


def write_points( points ):


    # Open directory
    
    directory = "lattice"
    
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    file = open(directory + "/points", 'w')


    
    # Number of points
    
    file.write( str(len(points)) )
    
    file.write( "\n" )



    # Write points

    for p in range( len(points) ):

        file.write( "%.0f " % points[p,0] )
        file.write( "%.0f " % points[p,1] )
        file.write( "%.0f\n" % points[p,2] )

    file.close()
