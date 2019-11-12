import salome

import sys

import numpy as np


def remove_neighbours( nb, pindex, n ):

    """
    New neighbour array. Remove non existing neighbours based on indices
    """

    newNeigh = -1 + np.zeros( (n, nb.shape[1]), dtype = np.int64 )

    q = len(nb[0])

    
    for i in range( len(pindex) ):
        
        if pindex[i] != -1:

            for vid in range(q):

                if nb[i,vid] != -1:

                    newNeigh[ pindex[i],vid ] = pindex[ int(nb[i,vid]) ]
    
    
    return newNeigh
