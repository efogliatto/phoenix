import sys

import os

import numpy

import SMESH, SALOMEDS

import salome

    
def write_vtk_cells( Mesh, lbmodel ):

    directory = "lattice"
    
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    file = open(directory + "/vtkCells", 'w')


    # Elements from Mesh

    elements = []

    if lbmodel.D() == 2:    
        elements = Mesh.GetElementsByType(SMESH.FACE)

    else:
        elements = Mesh.GetElementsByType(SMESH.VOLUME)
        

    file.write( str( len(elements) )  )
    file.write( " {}\n".format(len( Mesh.GetElemNodes( elements[0] ) )) )
    

    # Move over elements and write
    
    for el in elements:
  
        element_nodes = Mesh.GetElemNodes( el )


        # Write with VTK ordering

        if lbmodel.D() == 2:            

            file.write( "%.0f " % (element_nodes[0] - 1) )
            file.write( "%.0f " % (element_nodes[1] - 1) )
            file.write( "%.0f " % (element_nodes[3] - 1) )
            file.write( "%.0f " % (element_nodes[2] - 1) )            


        else:

            file.write( "%.0f " % (element_nodes[0] - 1) )
            file.write( "%.0f " % (element_nodes[3] - 1) )
            file.write( "%.0f " % (element_nodes[1] - 1) )
            file.write( "%.0f " % (element_nodes[2] - 1) )
            file.write( "%.0f " % (element_nodes[4] - 1) )
            file.write( "%.0f " % (element_nodes[7] - 1) )
            file.write( "%.0f " % (element_nodes[5] - 1) )
            file.write( "%.0f " % (element_nodes[6] - 1) )                        

            # for node in element_nodes:

            #     file.write( "%.0f " % (node-1) )
            
            
        file.write( "\n" )


    file.close()
