import sys

import os

import salome

import numpy

from salome.geom import geomtools
 
    
def write_vtk_cells( geompy, cells ):

    directory = "lattice"
    
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    file = open(directory + "/vtkCells", 'w')

    file.write( str( len(cells) )  )
    file.write( " 4\n" )
    

    for cell in cells:
        
        for i in cell:

            file.write( "%.0f " % i )

            
        file.write( "\n" )


    file.close()
