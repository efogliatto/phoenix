import salome

import sys

import numpy as np

from salome.geom import geomtools


def lattice_mesh_points(geompy, shape):
    
    '''
    This function returns lattice points inside the geometry bounding box
    
    No shape checking!

    '''

    bbox = geompy.BoundingBox( shape )

    nx = int(  bbox[1] - bbox[0]  ) + 1
    ny = int(  bbox[3] - bbox[2]  ) + 1
    

    pointList = np.zeros( (nx*ny, 3) )


    for j in range(ny):

        for i in range(nx):

            pointList[i+j*nx,0] = int( bbox[0] ) + i
            pointList[i+j*nx,1] = int( bbox[2] ) + j


        
    return pointList
