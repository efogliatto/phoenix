import salome

import sys

import numpy as np


def remove_points( points, active_points ):

    """
    New point array. Remove points based on active locations
    """

    # Total number of new points

    newp = np.sum(active_points)
    
    newPoints = np.zeros( (newp,3), dtype = np.int64)


    # Mapping indices: old to new indexing

    oldToNew = -1 + np.zeros( (len(points),1), dtype = np.int64 )

    id = 0

    for i,p in enumerate(points):

        if active_points[i] == 1:

          oldToNew[i] = id

          id = id + 1

          

    # Point array slicing

    id = 0
    
    for i in range( len(oldToNew) ):
        
        if oldToNew[i] != -1:
            
            newPoints[id,0] = points[i,0]
            newPoints[id,1] = points[i,1]
            newPoints[id,2] = points[i,2]            

            id = id + 1

            
    
    
    return newPoints, oldToNew
