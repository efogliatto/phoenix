#------------------ Generacion de bibliotecas ------------------#
 
# lattice-boltzmann equations library
add_library(tequations SHARED 
  TEquation.C
  simplifiedTEq.C
  )
