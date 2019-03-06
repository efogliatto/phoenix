import salome

import sys

import numpy as np

from salome.geom import geomtools
    

def vtk_cells(geompy, nb):

    vtk = []

    rev =  [0, 3, 4, 1, 2, 7, 8, 5, 6]

    
    for id in range( len(nb) ):

        cell = []
        
        cell.append( id )

        
        for k in [3,7,4]:

            aux = nb[id][k]
            
            if aux != -1:

                cell.append(aux)


                
        if len(cell) == 4:

            a = cell[3]

            cell[3] = cell[2]

            cell[2] = a

            vtk.append( cell )


            
    return vtk
