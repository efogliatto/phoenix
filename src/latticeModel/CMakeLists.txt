#------------------ Generacion de bibliotecas ------------------#
 
# latticeModel library
add_library(latticemodel SHARED 
  latticeModel.C
  D2Q9.C
  D3Q15.C
  latticeModelCreator.C
  )
