#------------------ Generacion de bibliotecas ------------------#
 
# lattice-boltzmann equations library
add_library(lbequations SHARED 
  lbEquation.C
  pseudopotential/pseudoPotEquation.C
  pseudopotential/LiMRTEq.C
  pseudopotential/XuMRTEq.C
  pseudopotential/PPEquation.C
  pseudopotential/pseudoPotEqHandler.C
  pseudopotential/boundary/ppBndCond.C
  pseudopotential/boundary/ppWetNodeBnd.C
  pseudopotential/boundary/ppNEBB.C
  pseudopotential/boundary/ppNEExt.C
  pseudopotential/boundary/ppFixedU.C
  pseudopotential/boundary/ppPeriodic.C
  pseudopotential/boundary/ppOutflow.C
  pseudopotential/boundary/ppBndCreator.C
  energy/energyEquation.C
  energy/myMRTEq.C
  energy/LiEnergyMRTEq.C
  energy/EEquation.C
  energy/energyEqHandler.C
  energy/boundary/energyBndCond.C
  energy/boundary/energyFixedT.C
  energy/boundary/energyFixedTSpot.C
  energy/boundary/energyFixedCosT.C
  energy/boundary/energyFixedGradT.C
  energy/boundary/energyPeriodic.C
  energy/boundary/energyNormalHeatFlux.C
  energy/boundary/energyOutflow.C    
  energy/boundary/energyBndCreator.C
  )
