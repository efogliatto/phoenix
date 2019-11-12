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
    

    
    for c in cells:

        if len(c) == 4:

            # Cell centre
            
            Vertex = geompy.MakeVertex( ( 0.5*points[ c[1],0 ] + 0.5*points[ c[0],0 ] ) , ( 0.5*points[ c[2],1 ] + 0.5*points[ c[1],1 ] ), 0 )

            Face = geompy.MakePlane(Vertex, OZ, 1)


            # Compute intersections
            
            Common = geompy.MakeCommonList([Face, shape], True)

            

            


        elif len(c) == 8:

            pass


    return cells
