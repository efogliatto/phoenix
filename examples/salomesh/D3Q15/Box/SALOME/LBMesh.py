# -*- coding: utf-8 -*-

###
### This file is generated automatically by SALOME v8.2.0 with dump python functionality
###

import sys

import salome

salome.salome_init()
theStudy = salome.myStudy

import salome_notebook
notebook = salome_notebook.NoteBook(theStudy)


# import own functions
sys.path.append('/users/fogliate/phoenix/src/tools/mesh')
import salomesh as sm



###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS
import numpy as np



# Box definition

dx = 2
dy = 2
dz = 2


geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )


Cavity = geompy.MakeBoxDXDYDZ(dx, dy, dz)

geompy.addToStudy( Cavity, 'Cavity' )




# List of faces

faceList = geompy.ExtractShapes(Cavity, geompy.ShapeType["FACE"], True)


# Groups

X0 = geompy.CreateGroup(Cavity, geompy.ShapeType["FACE"])
X1 = geompy.CreateGroup(Cavity, geompy.ShapeType["FACE"])
Y0 = geompy.CreateGroup(Cavity, geompy.ShapeType["FACE"])
Y1 = geompy.CreateGroup(Cavity, geompy.ShapeType["FACE"])
Z0 = geompy.CreateGroup(Cavity, geompy.ShapeType["FACE"])
Z1 = geompy.CreateGroup(Cavity, geompy.ShapeType["FACE"])


X0.SetName('X0')
X1.SetName('X1')
Y0.SetName('Y0')
Y1.SetName('Y1')
Z0.SetName('Z0')
Z1.SetName('Z1')


GroupsDict = {}
GroupsDict['X0'] = []
GroupsDict['X1'] = []
GroupsDict['Y0'] = []
GroupsDict['Y1'] = []
GroupsDict['Z0'] = []
GroupsDict['Z1'] = []




for face in faceList:

  com   = geompy.MakeCDG(face)
  
  coord = geompy.PointCoordinates(com)

  
  # X0
  
  if np.isclose(coord[0], 0, rtol=1e-05, atol=1e-08, equal_nan=False):

    geompy.UnionList(X0,[face])

    GroupsDict['X0'].append(face)


  # X1

  elif np.isclose(coord[0], dx, rtol=1e-05, atol=1e-08, equal_nan=False):

    geompy.UnionList(X1,[face])

    GroupsDict['X1'].append(face)

    
  # Y0

  elif np.isclose(coord[1], 0, rtol=1e-05, atol=1e-08, equal_nan=False): 

    geompy.UnionList(Y0,[face])

    GroupsDict['Y0'].append(face)
    
    
  # Y1

  elif np.isclose(coord[1], dy, rtol=1e-05, atol=1e-08, equal_nan=False):

    geompy.UnionList(Y1,[face])

    GroupsDict['Y1'].append(face)

  
  # Z0

  elif np.isclose(coord[2], 0, rtol=1e-05, atol=1e-08, equal_nan=False): 

    geompy.UnionList(Z0,[face])

    GroupsDict['Z0'].append(face)
    
    
  # Z1

  elif np.isclose(coord[2], dz, rtol=1e-05, atol=1e-08, equal_nan=False):

    geompy.UnionList(Z1,[face])

    GroupsDict['Z1'].append(face)

  

    
# List of groups

GroupsList = [X0,X1,Y0,Y1,Z0,Z1]

for group in GroupsList:

  geompy.addToStudyInFather( Cavity, group, group.GetName() )






# Integer bounding box

BBox = geompy.BoundingBox(Cavity, True)

Box = geompy.MakeBoxTwoPnt(  geompy.MakeVertex( np.ceil(BBox[0]), np.ceil(BBox[2]), np.ceil(BBox[4]) ),  geompy.MakeVertex( np.ceil(BBox[1]), np.ceil(BBox[3]), np.ceil(BBox[5])) )

geompy.addToStudy( Box, 'Bounding box' )



###
### SMESH component
###

import  SMESH, SALOMEDS

from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Local_Length = smesh.CreateHypothesis('LocalLength')

Local_Length.SetLength( 1 )

Local_Length.SetPrecision( 1e-07 )

Mesh = smesh.Mesh( Box )

Cartesian_3D = Mesh.BodyFitted()

Body_Fitting_Parameters = Cartesian_3D.SetGrid([ [ '1' ], [ 0, 1 ]],[ [ '1' ], [ 0, 1 ]],[ [ '1' ], [ 0, 1 ]],2,0)

Body_Fitting_Parameters.SetFixedPoint( SMESH.PointStruct ( 0, 0, 0 ), 1 )

Body_Fitting_Parameters.SetAxesDirs( SMESH.DirStruct( SMESH.PointStruct ( 1, 0, 0 )), SMESH.DirStruct( SMESH.PointStruct ( 0, 1, 0 )), SMESH.DirStruct( SMESH.PointStruct ( 0, 0, 1 )) )

isDone = Mesh.Compute()




# Filter elements lying on geometry
# This way can be used with multiple criterions

# criterion = smesh.GetCriterion(SMESH.VOLUME,SMESH.FT_LyingOnGeom,Sphere,SMESH.FT_LogicalNOT)
criterion = smesh.GetCriterion(SMESH.VOLUME,SMESH.FT_BelongToGeom,Cavity,SMESH.FT_LogicalNOT, Tolerance=0.7071)

filter = smesh.GetFilterFromCriteria([criterion])

isDone = Mesh.RemoveElements( Mesh.GetIdsFromFilter(filter) )

isDone = Mesh.RemoveOrphanNodes()

Mesh.RenumberNodes()

Mesh.RenumberElements()




# Create groups for boundary nodes


# Nodes connected to less than 8 elements are on boundary

bnd_criterion = smesh.GetCriterion(SMESH.NODE,SMESH.FT_NodeConnectivityNumber,SMESH.FT_LessThan,8)

bnd_filter = smesh.GetFilterFromCriteria([bnd_criterion])

bnd_ids = Mesh.GetIdsFromFilter(bnd_filter) 



# Create empty mesh groups based on geometry groups

Mesh_groups = {}

for group in GroupsList:

    Mesh_groups[group.GetName()] = Mesh.CreateEmptyGroup(SMESH.NODE, group.GetName())
    
    

# Look closest surface (group) and add to corresponding mesh group

for node in bnd_ids:

    node_xyz = Mesh.GetNodeXYZ( node )

    node_vertex = geompy.MakeVertex(node_xyz[0], node_xyz[1], node_xyz[2])

    dist = float(1000)

    closest_bnd = ""

    for group in GroupsList:

        distAux = geompy.MinDistance( group, node_vertex )

        # print('{} {} {}'.format(node, distAux, group.GetName()))        
        
        if distAux < dist:

            closest_bnd = group.GetName()

            dist = distAux
            
            
    Mesh_groups[ closest_bnd ].Add( [node] )

  





    




## Set names of Mesh objects
smesh.SetName(Cartesian_3D.GetAlgorithm(), 'Cartesian_3D')
smesh.SetName(Body_Fitting_Parameters, 'Body Fitting Parameters')
smesh.SetName(Local_Length, 'Local Length')
smesh.SetName(Mesh.GetMesh(), 'Mesh')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()






  
  


# ##############################
# #            MESH            #
# ##############################

# # Mesh creation

# mesh = sm.lbmesh(geompy, Cavity, lattice_model = "D3Q15")


# # Set boundaries from geometry groups

# mesh.setGroupsFromGeometry( [Z0,Z1,X0,X1,Y0,Y1] )


# # Set periodic boundaries

# mesh.setPeriodicBoundaries( [ ('X0', 'X1'), ('Y0', 'Y1') ] )


# # Set explicit periodic corner correction

# mesh.setCorners( [((0,0,0),(dx,0,0)), ((0,dy,0),(dx,dy,0)), ((0,0,0),(dx,0,dz)), ((0,dy,0),(dx,dy,dz))] )

# # mesh.setCorners( [((0,0,0),(dx,0,0)), ((0,dy,0),(dx,dy,0)), ((0,0,0),(dx,0,dz)), ((0,dy,0),(dx,dy,dz)), ((0,0,0),(0,dy,0)), ((dx,0,0),(dx,dy,0)), ((0,0,0),(0,dy,dz)), ((dx,0,0),(dx,dy,dz)), ((0,0,0),(dx,dy,0)), ((dx,0,0),(dx,dy,0)), ((0,0,dz),(dx,dy,dz)), ((dx,0,dz),(dx,dy,dz))] )


# # Mesh calculation and saving

# mesh.compute()

# mesh.export()
  




# # Mesh creation

# points, nb = sm.mesh(geompy, Cavity, "D2Q9")


# # VTKCells

# vtk = sm.vtk_cells(geompy, nb)


# # Boundaries detection

# bdWeights = {'X0': 1, 'X1': 1, 'Y0': 2, 'Y1': 2}

# bdDict = sm.lattice_boundaries_weights(geompy, Cavity, [X0,X1,Y0,Y1], points, nb, bdWeights)


# periodicBCs = [ ('X0', 'X1') ]

# # corners = [((x0,y0),(x3,y0)),
# #            ((x0,y1),(x3,y1))] 

# nb = sm.lattice_periodic_bcs(geompy, nb, bdDict, periodicBCs, points)





# # Write mesh info in LB format

# sm.write_points( geompy, points )

# sm.write_neighbours( geompy, nb )

# sm.write_boundaries( geompy, bdDict )

# sm.write_vtk_cells( geompy, vtk )




# if salome.sg.hasDesktop():
#   salome.sg.updateObjBrowser(True)
