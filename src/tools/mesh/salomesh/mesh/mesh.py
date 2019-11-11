import salome

import sys

import numpy as np

from salome.geom import geomtools

from .lattice_mesh_points import lattice_mesh_points

from .lattice_neighbours import lattice_neighbours

from .points_in_shape import points_in_shape

from .remove_neighbours import remove_neighbours


def mesh(geompy, shape, model):

    """
    Lattice mesh creation
    """


    # Point creation

    base_points = lattice_mesh_points(geompy, shape)


    
    # Neighbour creation

    base_nb = lattice_neighbours(geompy, shape, base_points, model)


    
    # Points inside geometry

    npid,points = points_in_shape(geompy, shape, base_points)


    
    # Non existing neighbours removal

    nb = remove_neighbours(geompy, npid, base_nb, len(points) )


    

    return points, nb
