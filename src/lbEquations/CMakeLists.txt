#------------------ Generacion de bibliotecas ------------------#
 
# lattice-boltzmann equations library
add_library(lbequations SHARED 
  lbEquation.C
  pseudopotential/pseudoPotEquation.C
  pseudopotential/LiMRTEq.C
  pseudopotential/PPEquation.C
  pseudopotential/pseudoPotEqHandler.C
  pseudopotential/boundary/ppBndCond.C
  pseudopotential/boundary/ppFixedU.C
  pseudopotential/boundary/ppPeriodic.C
  pseudopotential/boundary/ppBndCreator.C
  energy/energyEquation.C
  energy/myMRTEq.C
  energy/EEquation.C
  energy/energyEqHandler.C
  energy/boundary/energyBndCond.C
  energy/boundary/energyFixedT.C
  energy/boundary/energyFixedGradT.C
  energy/boundary/energyPeriodic.C
  energy/boundary/energyBndCreator.C
  
  )
