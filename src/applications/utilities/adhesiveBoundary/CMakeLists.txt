#------------------ Distributed adhesive coefficients for pseudopotential multiphase models ------------------#

add_executable(adhesiveBoundary
  "adhesiveBoundary.C"
  "spotSample/spotSample.C"
  "spotSample/uniformSpots.C"
  "spotSample/randomSpots.C"  
  "spotSample/spotSampleCreator.C"
  "spotRadius/spotRadius.C"
  "spotRadius/fixedRadiusSpot.C"
  "spotRadius/normalDistRadiusSpot.C"
  "spotRadius/spotRadiusCreator.C"
  )

target_link_libraries(adhesiveBoundary ${PROJECT_LINK_LIBS} ${MPI_LIBRARIES} ${PROJECT_LINK_LIBS})
