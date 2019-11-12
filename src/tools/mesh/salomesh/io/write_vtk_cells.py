import sys

import os

import numpy

    
def write_vtk_cells( cells ):

    directory = "lattice"
    
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    file = open(directory + "/vtkCells", 'w')

    file.write( str( len(cells) )  )
    file.write( " {}\n".format(len(cells[0])) )
    

    for cell in cells:
        
        for i in cell:

            file.write( "%.0f " % i )

            
        file.write( "\n" )


    file.close()
