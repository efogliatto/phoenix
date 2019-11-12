import salome

import sys

import numpy as np

from salome.geom import geomtools

from .lattice_mesh_points import lattice_mesh_points

from .lattice_neighbours import lattice_neighbours

from .points_in_shape import points_in_shape

from .remove_neighbours import remove_neighbours

from .vtk_cells import vtk_cells

from ..DdQq.DdQq import DdQq



class lbmesh:

    """
    Lattice boltzmann mesh class
    """

    def __init__(self, geompy, shape, lattice_model="D2Q9"):


        # Basic elements
        
        self.__lmodel = DdQq(lattice_model)

        self.__geompy = geompy

        self.__shape = shape

        pass



    def basicMesh(self):

        """
        Compute basic mesh 
        """

        # Point creation

        base_points, grid = lattice_mesh_points(self.__geompy, self.__shape, self.__lmodel.D())


    
        # Neighbour creation

        base_nb = lattice_neighbours(base_points, grid, self.__lmodel)



        # # VTKCells

        # base_vtk = vtk_cells(grid, self.__lmodel)

        
        
    
        # # Points inside geometry

        # npid,points = points_in_shape(geompy, shape, base_points)


    
        # # Non existing neighbours removal

        # nb = remove_neighbours(geompy, npid, base_nb, len(points) )


        pass



    def compute(self):
        """
        Compute lattice mesh
        """

        # Compute basic mesh
        self.basicMesh()
