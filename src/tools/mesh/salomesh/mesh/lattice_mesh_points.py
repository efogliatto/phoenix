import salome

import sys

import numpy as np

from salome.geom import geomtools


def lattice_mesh_points(geompy, shape, dim):
    
    '''
    This function returns lattice points inside the geometry bounding box
    
    No shape checking!

    '''


    # Compute bounding box
    
    bbox = geompy.BoundingBox( shape )

    nx = int(  bbox[1] - bbox[0]  ) + 1
    ny = int(  bbox[3] - bbox[2]  ) + 1
    nz = int(  bbox[5] - bbox[4]  ) + 1
    

    # Compute base mesh according to dimensions

    if dim == 2:   

        pointList = np.zeros( (nx*ny, 3) )

        for j in range(ny):

            for i in range(nx):

                pointList[i+j*nx,0] = int( bbox[0] ) + i
                pointList[i+j*nx,1] = int( bbox[2] ) + j



    elif dim == 3:

        pointList = np.zeros( (nx*ny*nz, 3) )

        for k in range(nz):
        
            for j in range(ny):

                for i in range(nx):

                    pointList[i + j*nx + k*nx*ny,0] = int( bbox[0] ) + i
                    pointList[i + j*nx + k*nx*ny,1] = int( bbox[2] ) + j
                    pointList[i + j*nx + k*nx*ny,2] = int( bbox[4] ) + k


    
    else:

        assert False, "Dimension error in lattice points construction"
    
        
    return pointList, (nx,ny,nz)
