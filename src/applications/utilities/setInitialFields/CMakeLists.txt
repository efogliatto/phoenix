#------------------ Initial fields setting ------------------#

add_executable(setInitialFields
  "setInitialFields.C"
  "scalarFields/initialScalarField.C"
  "scalarFields/boxBoundedScalarField.C"
  "scalarFields/bgSphereBoundedScalarField.C"  
  "scalarFields/uniformScalarField.C"
  "scalarFields/randomScalarField.C"
  "scalarFields/linearScalarField.C"
  "scalarFields/tgProfileScalarField.C"  
  )

target_link_libraries(setInitialFields ${PROJECT_LINK_LIBS} ${MPI_LIBRARIES} ${PROJECT_LINK_LIBS})
