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

  

    


for group in [X0,X1,Y0,Y1,Z0,Z1]:

  geompy.addToStudyInFather( Cavity, group, group.GetName() )





##############################
#       MESH CREATION        #
##############################

mesh = sm.lbmesh(geompy, Cavity)

mesh.compute()
  




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
