#------------------ Generacion de bibliotecas ------------------#
 
# LB Forces library
add_library(forces SHARED 
  EOS/EOS.C
  EOS/vanDerWaals.C
  EOS/CarnahanStarling.C
  EOS/EOSCreator.C
  interaction/interactionForce.C
  interaction/singleRangeIntForce.C
  interaction/singleRangeMixedIntForce.C
  interaction/intForce.C
  buoyant/buoyantForce.C
  buoyant/fixedDensityBForce.C
  buoyant/avgDensityBForce.C
  buoyant/bForces.C
  external/externalForce.C
  adhesive/adhesiveForce.C
  adhesive/noAds.C
  adhesive/liAdhesive.C
  adhesive/adsForceCreator.C
  pseudopotential/pseudoPotForce.C
  )
