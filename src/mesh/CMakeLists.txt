#------------------ Generacion de bibliotecas ------------------#
 
# Mesh library
add_library(mesh SHARED 
  latticeMesh.C
  localIndexer.C
  mpiInfo.C
  )
