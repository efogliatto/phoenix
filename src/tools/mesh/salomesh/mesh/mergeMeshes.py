import salome

import sys

import numpy as np

from salome.geom import geomtools

from ..DdQq.DdQq import DdQq

from ..io.write_points import write_points

from ..io.write_neighbours import write_neighbours

from ..io.write_vtk_cells import write_vtk_cells

from ..io.write_boundaries import write_boundaries

from .lbmesh import lbmesh


import SMESH, SALOMEDS

from salome.smesh import smeshBuilder



def mergeMeshes(geompy, mesh_a, mesh_b, computeNeighbours = True, name='Mesh'):
    
    """
    Merge two lbmesh    

    """

    merged = lbmesh(geompy, mesh_a.shape, mesh_a.lmodel.name())

    merged.smesh = mesh_a.smesh
    
    merged.Mesh = merged.smesh.Concatenate( [ mesh_a.Mesh.GetMesh(), mesh_b.Mesh.GetMesh() ], 1, 1, 1e-05, False )

    merged.smesh.SetName(merged.Mesh.GetMesh(), name)


    # Compute neighbours

    if computeNeighbours == True:

        merged.updateNeighbours()


    return merged

    
