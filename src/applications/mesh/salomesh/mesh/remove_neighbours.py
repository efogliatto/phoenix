import salome

import sys

import numpy as np

from salome.geom import geomtools


def remove_neighbours(geompy, pindex, nb, n ):

    """
    New neighbour array. Remove non existing neighbours based on indices
    """

    newNeigh = -1 + np.zeros( (n, nb.shape[1]), dtype = int )

    Q = nb.shape[1]


   
    for i, id in enumerate(pindex):
        
        if id != -1:

            for vid in range(Q):

                if nb[i][vid] != -1:

                    newNeigh[id][vid] = pindex[ nb[i][vid] ]
    
    
    return newNeigh
