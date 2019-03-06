import salome

import sys

import numpy as np

from salome.geom import geomtools


def lattice_neighbours(geompy, shape, points, model):
    
    '''
    This function returns neighbour indices
    '''

    nb = []

    bbox = geompy.BoundingBox( shape )
    
    if model is "D2Q9":


        
        # Resize neighbour

        nb = -1. + np.zeros( (len(points),9) )
            

        
        # Lattice velocities and reverse indices
        
        vel = np.array( [[0,0,0], [1,0,0], [0,1,0], [-1,0,0], [0,-1,0], [1,1,0], [-1,1,0], [-1,-1,0], [1,-1,0] ])

        rev =  [0, 3, 4, 1, 2, 7, 8, 5, 6]


        
        # Mesh. number of points
        
        nx = int(  bbox[1] - bbox[0]  ) + 1
        ny = int(  bbox[3] - bbox[2]  ) + 1
        nz = int(  bbox[5] - bbox[4]  ) + 1
        


        # Internal points first

        for j in range(1,ny-1):

            for i in range(1,nx-1):

                pointId = i + j*nx


                # Iterate on velocities

                for velId in range(9):

                    nb[pointId,velId] = pointId   +   vel[ rev[velId], 0]   +   vel[ rev[velId] , 1] * nx



                    
    # For boundary nodes, check neighbouring using distance to point

    for j in range(0,ny,ny-1):

        for i in range(nx):

            pointId = i+j*nx


    	    # Iterate on velocities

            for velId in range(9):
                
    		newId = pointId   +   vel[ rev[velId] ,0]   +   vel[ rev[velId] ,1] * nx
                

    		if  newId >= 0  and   newId <= nx*ny-1:

    		    if ( (  np.abs( points[pointId,0] - points[newId,0] ) <= 1  )   and   (  np.abs( points[pointId,1] - points[newId,1] ) <= 1  )  ):
	    
    			nb[pointId,velId] = newId



    for j in range(1,ny-1):

        for i in range(0,nx,nx-1):

    	    pointId = i + j*nx

            
     	    # Iterate on velocities
            
            for velId in range(9):

    		newId = pointId   +   vel[ rev[velId] ,0]   +   vel[ rev[velId] ,1] * nx
                

    		if  newId >= 0  and   newId <= nx*ny-1:

    		    if ( (  np.abs( points[pointId,0] - points[newId,0] ) <= 1  )   and   (  np.abs( points[pointId,1] - points[newId,1] ) <= 1  )  ):
	    
    			nb[pointId,velId] = newId
	        
        
        
    return nb
