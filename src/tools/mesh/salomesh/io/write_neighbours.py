import sys

import os

import numpy


def write_neighbours( nb ):

    directory = "lattice"
    
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    file = open(directory + "/neighbours", 'w')


    for nbi in nb:
        
        for i in nbi:

            file.write( "%.0f " % i )


            
        file.write( "\n" )


    file.close()
