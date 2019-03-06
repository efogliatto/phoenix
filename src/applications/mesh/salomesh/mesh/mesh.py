import salome

import sys

import numpy as np

from salome.geom import geomtools

import lattice_mesh_points as lmp

import lattice_neighbours as ln

import points_in_shape as pish

import remove_neighbours as rn


def mesh(geompy, shape, model):

    """
    Lattice mesh creation
    """


    # Point creation

    base_points = lmp.lattice_mesh_points(geompy, shape)


    
    # Neighbour creation

    base_nb = ln.lattice_neighbours(geompy, shape, base_points, model)


    
    # Points inside geometry

    npid,points = pish.points_in_shape(geompy, shape, base_points)


    
    # Non existing neighbours removal

    nb = rn.remove_neighbours(geompy, npid, base_nb, len(points) )


    

    return points, nb
