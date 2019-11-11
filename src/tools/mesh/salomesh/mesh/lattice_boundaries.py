import salome

import sys

import numpy as np

from salome.geom import geomtools


def lattice_boundaries(geompy, shape, ftList, points, nb):


    # Initialize boundary dictionary
    
    bdDict = {}

    for ft in ftList:

        bdDict[ ft.GetName() ] = []
        
    count = 0

    
    # Move over points. If some neighbour is -1, look for closest edge/face

    for id in range( len(points) ):

        
        bd = False

        
        for k in range( len(nb[id]) ):

            if nb[id][k] == -1:

                bd = True



       # Look for closest boundary

        if bd is True:

            
            dist = 1000

            ft_id = 0
        
        
            for ft in range( len(ftList) ):

                distAux = geompy.MinDistance( ftList[ft], points[id] )

                if distAux < dist:

                    dist = distAux

                    ft_id = ft



            # Add point index to boundary dictionary

            bdDict[ ftList[ft_id].GetName() ].append(id)



            

    return bdDict
