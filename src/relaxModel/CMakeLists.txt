#------------------ Relaxation factor model ------------------#
 
add_library(relaxmodel SHARED 
  relaxModel.C
  uniformTau.C
  rhoPieceWiseLinear.C
  rhoPieceWise.C
  relaxModelCreator.C
  )
