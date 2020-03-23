#------------------ Initial fields setting ------------------#

add_executable(setInitialFieldss
  "setInitialFields.C"
  "scalarFields/initialScalarField.C"
  "scalarFields/boxBoundedScalarField.C"  
  )

target_link_libraries(setInitialFieldss ${PROJECT_LINK_LIBS} ${MPI_LIBRARIES} ${PROJECT_LINK_LIBS})
