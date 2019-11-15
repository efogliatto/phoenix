import sys

import os

import numpy

    
def write_vtk_cells( Mesh ):

    directory = "lattice"
    
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    file = open(directory + "/vtkCells", 'w')


    # Elements from Mesh

    elements = Mesh.GetElementsId()    
    

    file.write( str( len(elements) )  )
    file.write( " {}\n".format(len( Mesh.GetElemNodes( elements[0] ) )) )
    

    # Move over elements and write
    
    for el in elements:
  
        element_nodes = Mesh.GetElemNodes( el )

        for node in element_nodes:

            file.write( "%.0f " % node )
            
            
        file.write( "\n" )


    file.close()
