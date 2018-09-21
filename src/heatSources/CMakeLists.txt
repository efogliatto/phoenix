#------------------ Generacion de bibliotecas ------------------#
 
# LB Heat Sources library
add_library(heatsources SHARED 
  heatSource.C
  markusHaziHS.C
  hsCreator.C
  )
