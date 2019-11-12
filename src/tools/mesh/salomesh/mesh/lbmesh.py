import salome

import sys

import numpy as np

from salome.geom import geomtools

from .lattice_mesh_points import lattice_mesh_points

from .lattice_neighbours import lattice_neighbours

from .points_in_shape import points_in_shape

from .remove_neighbours import remove_neighbours

from .remove_points import remove_points

from .vtk_cells import vtk_cells

from .remove_cells import remove_cells

from .update_vtkCells import update_vtkCells

from .check_active_points import check_active_points

from ..DdQq.DdQq import DdQq

from ..io.write_points import write_points

from ..io.write_neighbours import write_neighbours

from ..io.write_vtk_cells import write_vtk_cells



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


        
        pass


    
    def setCellFraction(f):
        """
        Volume fraction used in cell removal
        """

        self.__fraction = f

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

        
        
    
        # # Points inside geometry

        # npid,points = points_in_shape(geompy, shape, base_points)


    

        pass




    def compute(self):
        """
        Compute lattice mesh
        """

        # Compute basic mesh

        self.basicMesh()


        # Remove extra cells
        # Falta remover vecinos en caras convexas

        self.__vtkCells = remove_cells(self.__geompy, self.__shape, self.__points, self.__vtkCells, self.__fraction)


        # Activate points

        active_points = check_active_points(self.__vtkCells, len(self.__points))

        
        # Remove non-existing points

        self.__points, oldToNew = remove_points(self.__points, active_points)


        # Remove non-existing neighbours

        self.__nb = remove_neighbours( self.__nb, oldToNew, len(self.__points) )


        # Update vtkCells with new indexing

        self.__vtkCells = update_vtkCells(self.__vtkCells, oldToNew)
        
        

        pass
    
    


    def export(self):
        """
        Export lattice mesh
        """
        
        write_points(self.__points)

        write_neighbours(self.__nb)

        write_vtk_cells(self.__vtkCells)

        pass
