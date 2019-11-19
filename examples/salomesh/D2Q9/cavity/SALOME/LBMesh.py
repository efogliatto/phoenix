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

x0 = 0
x1 = 2
x2 = 3
x3 = 5
y0 = 0
y1 = 2
y2 = 7


# Total number of replicates

nrep = 1


geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )


Vertex_0 = geompy.MakeVertex(x0, y1, 0)
Vertex_1 = geompy.MakeVertex(x1, y1, 0)
Vertex_2 = geompy.MakeVertex(x1, y0, 0)
Vertex_3 = geompy.MakeVertex(x2, y0, 0)
Vertex_4 = geompy.MakeVertex(x2, y1, 0)
Vertex_5 = geompy.MakeVertex(x3, y1, 0)
Vertex_6 = geompy.MakeVertex(x3, y2, 0)
Vertex_7 = geompy.MakeVertex(x0, y2, 0)



# Cavity lines
Line_0 = geompy.MakeLineTwoPnt(Vertex_0, Vertex_1)
Line_1 = geompy.MakeLineTwoPnt(Vertex_1, Vertex_2)
Line_2 = geompy.MakeLineTwoPnt(Vertex_2, Vertex_3)
Line_3 = geompy.MakeLineTwoPnt(Vertex_3, Vertex_4)
Line_4 = geompy.MakeLineTwoPnt(Vertex_4, Vertex_5)
Line_5 = geompy.MakeLineTwoPnt(Vertex_5, Vertex_6)
Line_6 = geompy.MakeLineTwoPnt(Vertex_6, Vertex_7)
Line_7 = geompy.MakeLineTwoPnt(Vertex_7, Vertex_0)


# Make face from lines
Wire_1 = geompy.MakeWire([Line_0, Line_1, Line_2, Line_3, Line_4, Line_5, Line_6, Line_7], 1e-07)
Cavity = geompy.MakeFaceWires([Wire_1], 1)
CavityBase = geompy.MakeFaceWires([Wire_1], 1)


# Copy cavity
for i in range(nrep):
  Cavity = geompy.MakeFuseList([Cavity, geompy.MakeTranslation(CavityBase, x3*i, 0, 0)], True, True)


geompy.addToStudy( Cavity, 'Cavity' )
# geompy.addToStudy( CavityBase, 'CavityBase' )



# List of edges

edgeList = geompy.ExtractShapes(Cavity, geompy.ShapeType["EDGE"], True)


# Groups

X0 = geompy.CreateGroup(Cavity, geompy.ShapeType["EDGE"])
X1 = geompy.CreateGroup(Cavity, geompy.ShapeType["EDGE"])
Y0 = geompy.CreateGroup(Cavity, geompy.ShapeType["EDGE"])
Y1 = geompy.CreateGroup(Cavity, geompy.ShapeType["EDGE"])


X0.SetName('X0')
X1.SetName('X1')
Y0.SetName('Y0')
Y1.SetName('Y1')


GroupsDict = {}
GroupsDict['X0'] = []
GroupsDict['X1'] = []
GroupsDict['Y0'] = []
GroupsDict['Y1'] = []




for edge in edgeList:

  com   = geompy.MakeCDG(edge)
  
  coord = geompy.PointCoordinates(com)

  
  # X0
  
  if np.isclose(coord[0], x0, rtol=1e-05, atol=1e-08, equal_nan=False):

    geompy.UnionList(X0,[edge])

    GroupsDict['X0'].append(edge)


  # X1

  elif np.isclose(coord[0], x3*nrep, rtol=1e-05, atol=1e-08, equal_nan=False):

    geompy.UnionList(X1,[edge])

    GroupsDict['X1'].append(edge)
    
    
  # Y1

  elif np.isclose(coord[1], y2, rtol=1e-05, atol=1e-08, equal_nan=False):

    geompy.UnionList(Y1,[edge])

    GroupsDict['Y1'].append(edge)

    
  # Y0

  else:    

    geompy.UnionList(Y0,[edge])

    GroupsDict['Y0'].append(edge)



GroupsList = [Y0,Y1,X0,X1]    

for group in GroupsList:

  geompy.addToStudyInFather( Cavity, group, group.GetName() )





##############################
#            MESH            #
##############################

# Mesh creation

mesh = sm.lbmesh(geompy, Cavity, maxDim=(5,7,0))

mesh.setTolerance( 1e-03 )


# Mesh calculation: castelation from cartesian grid

isDone = mesh.compute()


# Group creation from geometry

mesh.GroupsFromGeometry(GroupsList)


# Assign periodic boundaries

mesh.PeriodicBoundaries( [('X0','X1','X')] )


# Connect neighbours from specific nodes

mesh.ForcePeriodicPoints( [((x0,y1,0),(x3,y1,0)), ((x0,y2,0),(x3,y2,0))] )


# Export mesh in LB format

mesh.export()  
