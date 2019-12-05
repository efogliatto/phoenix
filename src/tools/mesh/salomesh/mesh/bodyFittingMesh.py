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


class bodyFittingMesh:

    """
    Lattice boltzmann mesh class. 

    Compute using body fitting algorithm
    """

    def __init__(self, geompy, shape, lattice_model="D2Q9", maxDim=(0,0,0)):


        # Basic elements
        
        self.lmodel = DdQq(lattice_model)

        self.shape = shape

        self.geompy = geompy

        self.Tolerance = 0.707107

        self.mg = {}
        
        

        # SMESH Hypotesis

        if self.lmodel.D() == 3:
        
            self.smesh = smeshBuilder.New()

            Local_Length = self.smesh.CreateHypothesis('LocalLength')

            Local_Length.SetLength( 1 )

            Local_Length.SetPrecision( 1e-07 )

            self.Mesh = self.smesh.Mesh( self.shape )

            Cartesian_3D = self.Mesh.BodyFitted()

            Body_Fitting_Parameters = Cartesian_3D.SetGrid([ [ '1' ], [ 0, 1 ]],[ [ '1' ], [ 0, 1 ]],[ [ '1' ], [ 0, 1 ]],2,0)

            Body_Fitting_Parameters.SetFixedPoint( SMESH.PointStruct ( 0, 0, 0 ), 1 )

            Body_Fitting_Parameters.SetAxesDirs( SMESH.DirStruct( SMESH.PointStruct ( 1, 0, 0 )), SMESH.DirStruct( SMESH.PointStruct ( 0, 1, 0 )), SMESH.DirStruct( SMESH.PointStruct ( 0, 0, 1 )) )


        else:

            self.smesh = smeshBuilder.New()

            self.Mesh = self.smesh.Mesh( self.shape )

            Regular_1D = self.Mesh.Segment()

            Local_Length_1 = Regular_1D.LocalLength(1,None,1e-07)

            Quadrangle_2D = self.Mesh.Quadrangle(algo=smeshBuilder.QUADRANGLE)
            

        
        pass




    def setTolerance(self,tol=1):
        """
        Set mesh filter tolerance
        """

        self.Tolerance = tol

        pass


    
    

    def compute(self, name='Mesh'):
        """
        Compute lattice mesh
        """

        # Compute basic mesh

        print('\nComputing basic mesh')
        
        isDone = self.Mesh.Compute()



        # Remove non tetra/hexa

        print('\nRemoving non orthogonal elements')        

        if self.lmodel.D() == 3:
            
            criterion = self.smesh.GetCriterion(SMESH.VOLUME,SMESH.FT_EntityType,SMESH.FT_Undefined,SMESH.Entity_Hexa,SMESH.FT_LogicalNOT)

            filter = self.smesh.GetFilterFromCriteria([criterion])

            isDone = self.Mesh.RemoveElements( self.Mesh.GetIdsFromFilter(filter) )

            isDone = self.Mesh.RemoveOrphanNodes()

            self.Mesh.RenumberNodes()        

            self.Mesh.RenumberElements()                    

        else:
            
            criterion = self.smesh.GetCriterion(SMESH.ALL,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,self.shape,SMESH.FT_LogicalNOT,SMESH.FT_Undefined,self.Tolerance)
            # criterion = self.smesh.GetCriterion(SMESH.VOLUME,SMESH.FT_EntityType,SMESH.FT_Undefined,SMESH.Entity_Tetra,SMESH.FT_LogicalNOT)
            
            filter = self.smesh.GetFilterFromCriteria([criterion])
            
            isDone = self.Mesh.RemoveElements( self.Mesh.GetIdsFromFilter(filter) )

            isDone = self.Mesh.RemoveOrphanNodes()

            nbRemoved = self.Mesh.RenumberNodes()
                
            self.Mesh.RenumberElements()        
                 



        # Set names of Mesh objects

        self.smesh.SetName(self.Mesh.GetMesh(), name)



        # Neighbour calculation

        # print('\nComputing neighbours')

        # self.neighbours = np.zeros( (self.Mesh.NbNodes(), self.lmodel.Q()), dtype=np.int64 )

        # self.neighbours.fill(-1)


        # # Fill points array

        # self.points = np.zeros( (self.Mesh.NbNodes(), 3), dtype=np.int64 )

        # for node in self.Mesh.GetNodesId():

        #   xyz = self.Mesh.GetNodeXYZ( node )

        #   for j in range(3):
  
        #     self.points[node-1,j] = np.rint( xyz[j] )



        # # Element inspection for connectivity

        # elements = self.Mesh.GetElementsId()


        # # Reverse velocity index

        # reverse = self.lmodel.reverse()


        # for el in elements:
  
        #   element_nodes = self.Mesh.GetElemNodes( el )

        #   dist = np.array([0,0,0], dtype=np.int64)
          
        #   for node in element_nodes:

        #     for nbnode in element_nodes:

        #       for j in range(3):

        #         dist[j] = np.rint( self.points[nbnode-1,j] - self.points[node-1,j] )


        #       # Corresponding velocity index

        #       vid = self.lmodel.vindex( dist )
              

        #       # Assign neighbour using reverse indexing

        #       if vid != -1:
        #           self.neighbours[node-1, reverse[vid]] = nbnode - 1

              

                
        return isDone
    
    


    

    def GroupsFromGeometry(self, geoGroups):
        """
        Create node groups from geometry groups
        """

        print('\nBoundary assignment')

        
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

          # node_xyz = self.Mesh.GetNodeXYZ( node )
          node_xyz = self.Mesh.GetNodeXYZ( node )
            
          # node_vertex = self.geompy.MakeVertex(node_xyz[0], node_xyz[1], node_xyz[2])
          node_vertex = self.geompy.MakeVertex(np.float64( node_xyz[0]), np.float64(node_xyz[1]), np.float64(node_xyz[2]) )          

          dist = float(1000)

          closest_bnd = ""

          for group in geoGroups:

            distAux = self.geompy.MinDistance( group, node_vertex )
        
            if distAux < dist:

              closest_bnd = group.GetName()

              dist = distAux
            
            
          Mesh_groups[ closest_bnd ].Add( [node] )

          self.mg[closest_bnd].append(node)





    def PeriodicBoundaries(self, periodic):
        """
        Change neighbours for periodic boundaries
        
        Arguments
        periodic: list of tuples ('bd1','bd2','dir'), where dir can be X,Y,Z
        """


        for pair in periodic:


            # Node indices for each boundary
            
            bd1 = self.mg[ pair[0] ]

            bd2 = self.mg[ pair[1] ]

            dir = pair[2]


            for nd1 in bd1:

                # pos_1 = np.rint(  self.Mesh.GetNodeXYZ( nd1 )  )
                pos_1 = self.points[ nd1-1 ]
                
                for nd2 in bd2:

                    # pos_2 = np.rint(  self.Mesh.GetNodeXYZ( nd2 )  )
                    pos_2 = self.points[ nd2-1  ]


                    if dir == 'X':

                        if (pos_1[1] == pos_2[1]):

                            if (pos_1[2] == pos_2[2]):

                                
                                for k in range( self.lmodel.Q() ):

                                    if self.neighbours[nd2-1,k] == -1:
                                        self.neighbours[nd2-1,k] = self.neighbours[nd1-1,k]

                                    if self.neighbours[nd1-1,k] == -1:
                                        self.neighbours[nd1-1,k] = self.neighbours[nd2-1,k]
                            

                    elif dir == 'Y':

                        if (pos_1[0] == pos_2[0]):

                            if (pos_1[2] == pos_2[2]):

                                for k in range( self.lmodel.Q() ):

                                    if self.neighbours[nd2-1,k] == -1:
                                        self.neighbours[nd2-1,k] = self.neighbours[nd1-1,k]

                                    if self.neighbours[nd1-1,k] == -1:
                                        self.neighbours[nd1-1,k] = self.neighbours[nd2-1,k]                                
                            

                    elif dir == 'Z':

                        if (pos_1[0] == pos_2[0]):

                            if (pos_1[1] == pos_2[1]):

                                
                                for k in range( self.lmodel.Q() ):

                                    if self.neighbours[nd2-1,k] == -1:
                                        self.neighbours[nd2-1,k] = self.neighbours[nd1-1,k]

                                    if self.neighbours[nd1-1,k] == -1:
                                        self.neighbours[nd1-1,k] = self.neighbours[nd2-1,k]                                
                    


        
        pass


          

    def ForcePeriodicPoints(self, pairsList):
        """
        Force neighbouring between pairs of points
        """

        for pair in pairsList:

            pt_1 = pair[0]

            pt_2 = pair[1]


            # Find node indices using closest node

            nd1 = self.Mesh.FindNodeClosestTo( pt_1[0], pt_1[1], pt_1[2] )

            nd2 = self.Mesh.FindNodeClosestTo( pt_2[0], pt_2[1], pt_2[2] )


            for k in range( self.lmodel.Q() ):

                if self.neighbours[nd2-1,k] == -1:
                    self.neighbours[nd2-1,k] = self.neighbours[nd1-1,k]

                if self.neighbours[nd1-1,k] == -1:
                    self.neighbours[nd1-1,k] = self.neighbours[nd2-1,k]
            
        
        pass
    

    
    def export(self):
        """
        Export lattice mesh
        """

        # Clear directory
        
        os.system( 'rm -rf lattice' )

        
        write_points(self.Mesh)

        # write_neighbours(self.neighbours)

        write_vtk_cells(self.Mesh, self.lmodel)

        write_boundaries(self.mg)

        pass
