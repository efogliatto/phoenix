import salome

import sys

import numpy as np

from salome.geom import geomtools


def lattice_neighbours(points, grid, lmodel):
    
    '''
    This function returns neighbour indices
    '''
    
        
    # Lattice velocities and reverse indices
        
    vel = lmodel.velocities()

    rev = lmodel.reverse()

    q = len(rev)

        
    # Mesh. number of points
        
    nx = int( grid[0] )
    ny = int( grid[1] )
    nz = int( grid[2] )


    # Resize neighbour

    nb = -1. + np.zeros( (len(points), q) )
    

    
    
    if lmodel.D() == 2:                


        for j in range(ny):

            for i in range(nx):

                pointId = i + j*nx



                # Boundary nodes

                if( (i == 0)  or  (i == nx-1)   or   (j == 0)  or  (j == ny-1) ):


                    # Iterate on velocities

                    for velId in range(q):
                
                        newId = pointId   +   vel[ rev[velId] ,0]   +   vel[ rev[velId] ,1] * nx
                    
                        if  newId >= 0  and   newId <= nx*ny-1:

                            if ( (  np.abs( points[pointId,0] - points[newId,0] ) <= 1  )   and   (  np.abs( points[pointId,1] - points[newId,1] ) <= 1  ) ):

                                nb[pointId,velId] = newId



                # Internal nodes

                else:

                    for velId in range(q):

                        nb[pointId,velId] = pointId   +   vel[ rev[velId], 0]   +   vel[ rev[velId] , 1] * nx
                



        

        # # Internal points first

        # for j in range(1,ny-1):

        #     for i in range(1,nx-1):

        #         pointId = i + j*nx


        #         # Iterate on velocities

        #         for velId in range(q):

        #             nb[pointId,velId] = pointId   +   vel[ rev[velId], 0]   +   vel[ rev[velId] , 1] * nx



                    
        # # For boundary nodes, check neighbouring using distance to point

        # for j in range(0,ny,ny-1):

        #     for i in range(nx):

        #         pointId = i+j*nx


    	#         # Iterate on velocities

        #         for velId in range(q):
                
        #             newId = pointId   +   vel[ rev[velId] ,0]   +   vel[ rev[velId] ,1] * nx
                

        #             if newId >= 0  and   newId <= nx*ny-1:

        #                 if ( (  np.abs( points[pointId,0] - points[newId,0] ) <= 1  )   and   (  np.abs( points[pointId,1] - points[newId,1] ) <= 1  )  ):
	    
        #                     nb[pointId,velId] = newId

                        


        # for j in range(1,ny-1):

        #     for i in range(0,nx,nx-1):

        #         pointId = i + j*nx

            
     	#         # Iterate on velocities
            
        #         for velId in range(q):

        #             newId = pointId   +   vel[ rev[velId] ,0]   +   vel[ rev[velId] ,1] * nx
                

        #             if  newId >= 0  and   newId <= nx*ny-1:

        #                 if ( (  np.abs( points[pointId,0] - points[newId,0] ) <= 1  )   and   (  np.abs( points[pointId,1] - points[newId,1] ) <= 1  )  ):
	    
        #                     nb[pointId,velId] = newId
            






    elif lmodel.D() == 3:
             

        for k in range(nz):

            for j in range(ny):

                for i in range(nx):

                    pointId = i + j*nx + k*nx*ny


                    # Boundary nodes

                    if( (i == 0)  or  (i == nx-1)   or   (j == 0)  or  (j == ny-1)   or   (k == 0)  or  (k == nz-1) ):

                        
                        # Iterate on velocities

                        for velId in range(q):
                
                            newId = pointId   +   vel[ rev[velId] ,0]   +   vel[ rev[velId] ,1] * nx   +   vel[ rev[velId], 2] * nx * ny
                    
                            if  newId >= 0  and   newId <= nx*ny*nz-1:

                                if ( (  np.abs( points[pointId,0] - points[newId,0] ) <= 1  )   and   (  np.abs( points[pointId,1] - points[newId,1] ) <= 1  )      and   (  np.abs( points[pointId,2] - points[newId,2] ) <= 1  )  ):

                                    nb[pointId,velId] = newId



                    # Internal nodes

                    else:

                        for velId in range(q):

                            nb[pointId,velId] = pointId   +   vel[ rev[velId], 0]   +   vel[ rev[velId] , 1] * nx   +   vel[ rev[velId] , 2] * nx*ny
                                    

        
    
                        
        
    return nb
