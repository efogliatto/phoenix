import salome

import sys

import numpy as np

from salome.geom import geomtools

from ..DdQq.DdQq import DdQq

from ..io.write_points import write_points

from ..io.write_neighbours import write_neighbours

from ..io.write_vtk_cells import write_vtk_cells

from ..io.write_boundaries import write_boundaries


import SMESH, SALOMEDS

from salome.smesh import smeshBuilder


class lbmesh:

    """
    Lattice boltzmann mesh class
    """

    def __init__(self, geompy, shape, lattice_model="D2Q9"):


        # Basic elements
        
        self.lmodel = DdQq(lattice_model)

        self.shape = shape

        self.geompy = geompy
        


        # Integer bounding box

        BBox = geompy.BoundingBox(shape, True)

        Box = geompy.MakeBoxTwoPnt(  geompy.MakeVertex( np.ceil(BBox[0]), np.ceil(BBox[2]), np.ceil(BBox[4]) ),  geompy.MakeVertex( np.ceil(BBox[1]), np.ceil(BBox[3]), np.ceil(BBox[5])) )

        geompy.addToStudy( Box, 'Bounding box' )

        

        # SMESH Hypotesis
        
        self.smesh = smeshBuilder.New()

        Local_Length = self.smesh.CreateHypothesis('LocalLength')

        Local_Length.SetLength( 1 )

        Local_Length.SetPrecision( 1e-07 )

        self.Mesh = self.smesh.Mesh( Box )

        Cartesian_3D = self.Mesh.BodyFitted()

        Body_Fitting_Parameters = Cartesian_3D.SetGrid([ [ '1' ], [ 0, 1 ]],[ [ '1' ], [ 0, 1 ]],[ [ '1' ], [ 0, 1 ]],2,0)

        Body_Fitting_Parameters.SetFixedPoint( SMESH.PointStruct ( 0, 0, 0 ), 1 )

        Body_Fitting_Parameters.SetAxesDirs( SMESH.DirStruct( SMESH.PointStruct ( 1, 0, 0 )), SMESH.DirStruct( SMESH.PointStruct ( 0, 1, 0 )), SMESH.DirStruct( SMESH.PointStruct ( 0, 0, 1 )) )
        

        
        pass


    

    def compute(self):
        """
        Compute lattice mesh
        """

        # Compute basic mesh

        isDone = self.Mesh.Compute()


        # Filter elements lying on geometry
        # This way can be used with multiple criterions

        criterion = self.smesh.GetCriterion(SMESH.VOLUME,SMESH.FT_BelongToGeom,self.shape,SMESH.FT_LogicalNOT, Tolerance=0.7071)

        filter = self.smesh.GetFilterFromCriteria([criterion])

        isDone = self.Mesh.RemoveElements( self.Mesh.GetIdsFromFilter(filter) )

        isDone = self.Mesh.RemoveOrphanNodes()

        self.Mesh.RenumberNodes()

        self.Mesh.RenumberElements()        



        # ## Set names of Mesh objects
        # smesh.SetName(Cartesian_3D.GetAlgorithm(), 'Cartesian_3D')
        # smesh.SetName(Body_Fitting_Parameters, 'Body Fitting Parameters')
        # smesh.SetName(Local_Length, 'Local Length')
        self.smesh.SetName(self.Mesh.GetMesh(), 'Mesh')



        # Neighbour calculation

        self.neighbours = np.zeros( (self.Mesh.NbNodes(), self.lmodel.Q()), dtype=np.int64 )

        self.neighbours.fill(-1)


        # Fill points array

        points = np.zeros( (self.Mesh.NbNodes(), 3), dtype=np.int64 )

        for node in self.Mesh.GetNodesId():

          xyz = self.Mesh.GetNodeXYZ( node )

          for j in range(3):
  
            points[node-1,j] = int( xyz[j] )



        # Element inspection for connectivity

        elements = self.Mesh.GetElementsId()


        # Reverse velocity index

        reverse = self.lmodel.reverse()


        for el in elements:
  
          element_nodes = self.Mesh.GetElemNodes( el )

          dist = np.array([0,0,0], dtype=np.int64)

          for node in element_nodes:

            for nbnode in element_nodes:

              for j in range(3):

                dist[j] = points[nbnode-1,j] - points[node-1,j]


              # Corresponding velocity index

              vid = self.lmodel.vindex( dist )


              # Assign neighbour using reverse indexing

              self.neighbours[node-1, reverse[vid]] = nbnode - 1
              

                
        return isDone
    
    


    

    def GroupsFromGeometry(self, geoGroups):
        """
        Create node groups from geometry groups
        """

        # Create groups for boundary nodes
        # Nodes connected to less than 8 (3D) or 4 (2D) elements are on boundary


        ncon = 8

        if self.lmodel.D() == 2:
            ncon = 4

            
        
        bnd_criterion = self.smesh.GetCriterion(SMESH.NODE,SMESH.FT_NodeConnectivityNumber,SMESH.FT_LessThan,ncon)

        bnd_filter = self.smesh.GetFilterFromCriteria([bnd_criterion])

        bnd_ids = self.Mesh.GetIdsFromFilter(bnd_filter) 



        # Create empty mesh groups based on geometry groups

        Mesh_groups = {}

        self.mg = {}

        for group in geoGroups:

          Mesh_groups[group.GetName()] = self.Mesh.CreateEmptyGroup(SMESH.NODE, group.GetName())

          self.mg[group.GetName()] = []
    
    

        # Look closest surface (group) and add to corresponding mesh group

        for node in bnd_ids:

          node_xyz = self.Mesh.GetNodeXYZ( node )

          node_vertex = self.geompy.MakeVertex(node_xyz[0], node_xyz[1], node_xyz[2])

          dist = float(1000)

          closest_bnd = ""

          for group in geoGroups:

            distAux = self.geompy.MinDistance( group, node_vertex )
        
            if distAux < dist:

              closest_bnd = group.GetName()

              dist = distAux
            
            
          Mesh_groups[ closest_bnd ].Add( [node] )

          self.mg[closest_bnd].append(node)



        
    

    def export(self):
        """
        Export lattice mesh
        """
        
        write_points(self.Mesh)

        write_neighbours(self.neighbours)

        write_vtk_cells(self.Mesh)

        write_boundaries(self.mg)

        pass
