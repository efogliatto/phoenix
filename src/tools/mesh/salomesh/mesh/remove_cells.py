import salome

import sys

import numpy as np

from salome.geom import geomtools


def remove_cells(geompy, shape, points, cells, fraction = 0.5):

    """
    Check cell intersection with shape
    """

    OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)


    # Cell removal array

    rm = np.zeros( (len(cells),1), dtype=np.int32 )
    

    if len(cells[0]) == 4:

        
        for i,c in enumerate(cells):
        

            # Cell centre
            
            Vertex = geompy.MakeVertex( ( 0.5*points[ c[1],0 ] + 0.5*points[ c[0],0 ] ) , ( 0.5*points[ c[2],1 ] + 0.5*points[ c[1],1 ] ), 0 )

            Face = geompy.MakePlane(Vertex, OZ, 1)


            # Compute intersections
            
            Common = geompy.MakeCommonList([Face, shape], True)

            prop = geompy.BasicProperties(Common)

            if prop[1] >= fraction:
                rm[i] = 1

            


    elif len(cells[0]) == 8:

        
        for i,c in enumerate(cells):

            
            # Cell vertices
            
            Vertex_0 = geompy.MakeVertex( points[ c[0],0 ], points[ c[0],1 ], points[ c[0],2 ] )

            Vertex_1 = geompy.MakeVertex( points[ c[7],0 ], points[ c[7],1 ], points[ c[7],2 ] )

            Box = geompy.MakeBoxTwoPnt(Vertex_0, Vertex_1)


            # Compute intersections
            
            Common = geompy.MakeCommonList([Box, shape], True)

            prop = geompy.BasicProperties(Common)

            if prop[2] >= fraction:
                rm[i] = 1            



    # Count number of active cells

    active = np.sum(rm)

    newCells = np.zeros( (active, len(cells[0])), dtype=np.int64 )


    # Copy only active cells

    nc = 0
    
    for i,c in enumerate(cells):

        if rm[i] == 1:
                
            newCells[nc] = c

        nc = nc + 1

    

    return newCells
