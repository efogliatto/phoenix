import salome

import sys

import numpy as np

from salome.geom import geomtools


def points_in_shape(geompy, shape, points):

    """
    Check which points are inside shape, and change index
    """

    newPid = -1 + np.zeros( len(points), dtype = int )

    newPoints = []

    id = 0


    for i in range( len(points) ):

        pt = geompy.MakeVertex(points[i,0], points[i,1], points[i,2])

        if np.isclose( geompy.MinDistance(shape, pt), 0 ) == True:

            newPid[i] = id
            
            id = id + 1

            newPoints.append(pt)
    
    

    return newPid, newPoints
