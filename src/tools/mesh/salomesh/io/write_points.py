import sys

import os

import salome

import numpy

# from salome.geom import geomtools


def write_points( Mesh ):


    # Open directory
    
    directory = "lattice"
    
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    file = open(directory + "/points", 'w')


    
    # Number of points
    
    file.write( '{}'.format( Mesh.NbNodes() ) )
    
    file.write( "\n" )


    # Write points
    
    for node in Mesh.GetNodesId():

      xyz = Mesh.GetNodeXYZ( node )

      file.write( "%.0f " % xyz[0] )

      file.write( "%.0f " % xyz[1] )

      file.write( "%.0f\n" % xyz[2] )


    file.close()
