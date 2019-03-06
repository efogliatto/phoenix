import salome

import sys

import numpy as np

from salome.geom import geomtools


def lattice_boundaries_weights(geompy, shape, GroupList, points, nb, bdWeights):


            

    # Initialize boundary dictionary
    
    bdDict = {}

    for ft in GroupList:

        bdDict[ ft.GetName() ] = []
        
    count = 0


    
    # Move over points. If some neighbour is -1, look for closest edge/face

    for id in range( len(points) ):

        
        bd = False

        
        for k in range( len(nb[id]) ):

            if nb[id][k] == -1:

                bd = True


       # Look for closest boundary.
       # Preference for periodic boundaries

        if bd == True:

            
            dist = []
                    
            for ft in GroupList:

                distAux = geompy.MinDistance( ft, points[id] )

                dist.append( (ft.GetName(), distAux) )



                
            # Sort according to distance

            sortedDist = sorted( dist, key=lambda ftDist: ftDist[1] )



            
            # First check if two features have same distance. Then assign priority for boundary elements

            if np.isclose( sortedDist[0][1], sortedDist[1][1], rtol=1e-05, atol=1e-08, equal_nan=False) == True:

                if bdWeights[ sortedDist[0][0] ] < bdWeights[ sortedDist[1][0] ]:
                    
                    c = sortedDist[0]

                    sortedDist[0] = sortedDist[1]

                    sortedDist[1] = c




            bdDict[ sortedDist[0][0] ].append(id)



            

    return bdDict
