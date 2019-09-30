#------------------ Distributed adhesive coefficients for pseudopotential multiphase models ------------------#

add_executable(adhesiveBoundary
  "adhesiveBoundary.C"
  "spotSample.C"
  "uniformSpots.C"
  "randomSpots.C"  
  "spotSampleCreator.C" )

target_link_libraries(adhesiveBoundary ${PROJECT_LINK_LIBS} ${MPI_LIBRARIES} ${PROJECT_LINK_LIBS})
