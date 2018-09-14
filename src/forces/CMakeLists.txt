#------------------ Generacion de bibliotecas ------------------#
 
# LB Forces library
add_library(forces SHARED 
  EOS.C
  vanDerWaals.C
  EOSCreator.C
  interactionForce.C
  singleRangeIntForce.C
  intForce.C
  buoyantForce.C
  fixedDensityBForce.C
  bForces.C
  )
