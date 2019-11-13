import sys

import os

import numpy


def write_boundaries( bd ):

    directory = "lattice"
    
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    file = open(directory + "/boundary", 'w')


    # Number of points
    
    file.write( str(len(bd.keys())) )
    file.write( "\n" )


    for key, value in bd.items() :

        file.write( "\n" )
        file.write( key )
        file.write( "\n" )
        file.write( str(len(value)) )
        file.write( "\n" )


        for id in value:
            file.write("%d\n" % id)
        
    

    file.close()
