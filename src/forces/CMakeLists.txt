#------------------ Generacion de bibliotecas ------------------#
 
# LB Forces library
add_library(forces SHARED 
  EOS.C
  vanDerWaals.C
  CarnahanStarling.C
  EOSCreator.C
  interactionForce.C
  singleRangeIntForce.C
  singleRangeMixedIntForce.C
  intForce.C
  buoyantForce.C
  fixedDensityBForce.C
  avgDensityBForce.C
  bForces.C
  externalForce.C
  pseudoPotForce.C
  )
