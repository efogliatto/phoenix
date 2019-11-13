import salome

import sys

import numpy as np

from salome.geom import geomtools

import time

from .lattice_mesh_points import lattice_mesh_points

from .lattice_neighbours import lattice_neighbours

from .points_in_shape import points_in_shape

from .remove_neighbours import remove_neighbours

from .remove_points import remove_points

from .vtk_cells import vtk_cells

from .lattice_boundaries_weights import lattice_boundaries_weights

from .lattice_periodic_bcs import lattice_periodic_bcs

from .remove_cells import remove_cells

from .update_vtkCells import update_vtkCells

from .check_active_points import check_active_points

from ..DdQq.DdQq import DdQq

from ..io.write_points import write_points

from ..io.write_neighbours import write_neighbours

from ..io.write_vtk_cells import write_vtk_cells

from ..io.write_boundaries import write_boundaries



class lbmesh:

    """
    Lattice boltzmann mesh class
    """

    def __init__(self, geompy, shape, lattice_model="D2Q9"):


        # Basic elements
        
        self.__lmodel = DdQq(lattice_model)

        self.__geompy = geompy

        self.__shape = shape
       

        self.__fraction = np.float64(0.5)

        self.__boundaries = {}

        self.__periodic = []

        self.__corners = []

        
        pass


    
    def setCellFraction(self,f):
        """
        Volume fraction used in cell removal
        """

        self.__fraction = f

        pass



    def setGroupsFromGeometry(self,groups):
        """
        Boundary assignment using geometry groups
        """

        self.__gfg = groups

        pass


    def setPeriodicBoundaries(self, periodic):
        """
        Pairs of periodic boundaries
        """

        self.__periodic = periodic

        pass


    def setCorners(self, corners):
        """
        Explicit corner correction.
        """

        self.__corners = corners

        pass

    

    def basicMesh(self):

        """
        Compute basic mesh 
        """

        # Point creation

        self.__points, grid = lattice_mesh_points(self.__geompy, self.__shape, self.__lmodel.D())


    
        # Neighbour creation

        self.__nb = lattice_neighbours(self.__points, grid, self.__lmodel)



        # VTKCells

        self.__vtkCells = vtk_cells(grid, self.__lmodel)

    

        pass




    def compute(self):
        """
        Compute lattice mesh
        """

        # Compute basic mesh

        start = time.time()
        self.basicMesh()
        end = time.time()

        print('\nBasic mesh computed in {:.4f} seconds\n'.format(end-start))


        # Remove extra cells

        start = time.time()
        self.__vtkCells, self.__nb = remove_cells(self.__geompy, self.__shape, self.__points, self.__nb, self.__vtkCells, self.__fraction)
        end = time.time()
        print('\nExtra cells removed in {:.4f} seconds\n'.format(end-start))
        
        
        # # Activate points

        # active_points = check_active_points(self.__vtkCells, len(self.__points))

        
        # # Remove non-existing points

        # self.__points, oldToNew = remove_points(self.__points, active_points)


        # # Remove non-existing neighbours

        # self.__nb = remove_neighbours( self.__nb, oldToNew, len(self.__points) )


        # # Update vtkCells with new indexing

        # self.__vtkCells = update_vtkCells(self.__vtkCells, oldToNew)


        # # Assign boundaries

        # self.__boundaries = lattice_boundaries_weights(self.__geompy, self.__shape, self.__gfg, self.__points, self.__nb)


        # # Assign periodic boundaries

        # if self.__periodic:
        #     self.__nb = lattice_periodic_bcs(self.__nb, self.__boundaries, self.__periodic, self.__points, self.__corners)
        

        pass
    
    


    def export(self):
        """
        Export lattice mesh
        """
        
        write_points(self.__points)

        write_neighbours(self.__nb)

        write_vtk_cells(self.__vtkCells)

        write_boundaries(self.__boundaries)

        pass
