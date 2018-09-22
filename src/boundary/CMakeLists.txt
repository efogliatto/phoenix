#------------------ Generacion de bibliotecas ------------------#
 
# Boundary conditions library
add_library(boundary SHARED 
  pseudopotential/ppBndCond.C
  pseudopotential/ppFixedU.C
  pseudopotential/ppPeriodic.C
  pseudopotential/ppBndCreator.C
  )
