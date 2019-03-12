#------------------ Generacion de bibliotecas ------------------#
 
# LB Forces library
add_library(forces SHARED 
  EOS/EOS.C
  EOS/vanDerWaals.C
  EOS/CarnahanStarling.C
  EOS/PengRobinson.C
  EOS/piecewiseLinear.C
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
  adhesive/phiBasedMod.C
  adhesive/adsForceCreator.C
  pseudopotential/pseudoPotForce.C
  surfaceTension/surfaceTension.C
  surfaceTension/noSurfaceTension.C
  surfaceTension/stCreator.C
  surfaceTension/liSurfaceTension.C  
  )
