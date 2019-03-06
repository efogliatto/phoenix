import salome

import sys

import numpy as np

from salome.geom import geomtools


def lattice_periodic_bcs(geompy, nb, bdDict, periodicBCs, points, corners = []):


    newNb = np.copy(nb)
    

    # Move over periodic boundaries and reassign neighbours

    for bdPair in periodicBCs:

        
        for i in range(  len( bdDict[ bdPair[0] ] )  ):

            id = bdDict[ bdPair[0] ][i]

            nid = bdDict[ bdPair[1] ][i]
                

            for k in range(9):

                
                if newNb[id][k] == -1:

                    newNb[id][k] = newNb[nid][k]

                    
                if newNb[nid][k] == -1:

                    newNb[nid][k] = newNb[id][k]

                    


    if len(corners) != 0:


        # Detect corners id's
        
        pairId = []
        
        for pair in corners:

            ids = [0,0]
            
            for i,p in enumerate(pair):

                
                # Move over points

                for pt in range(len(points)):

                    coord = geompy.PointCoordinates(points[pt])
                    
                    if coord[0] == p[0]  and  coord[1] == p[1]:

                        ids[i] = pt

                        
            pairId.append( ids )


            


        for pair in pairId:

            for k in range(9):

                if newNb[pair[0]][k] == -1:

                    newNb[pair[0]][k] = newNb[pair[1]][k]
            

                if newNb[pair[1]][k] == -1:

                    newNb[pair[1]][k] = newNb[pair[0]][k]                    


                    
    return newNb
