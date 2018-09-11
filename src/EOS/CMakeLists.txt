#------------------ Generacion de bibliotecas ------------------#
 
# EOS library
add_library(eos SHARED 
  EOS.C
  vanDerWaals.C
  EOSCreator.C
  )
