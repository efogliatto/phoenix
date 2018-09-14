#------------------ Generacion de bibliotecas ------------------#
 
# lattice-boltzmann equations library
add_library(lbequations SHARED 
  lbEquation.C
  pseudopotential/pseudoPotEquation.C
  pseudopotential/LiMRTEq.C
  pseudopotential/PPEquation.C
  )
