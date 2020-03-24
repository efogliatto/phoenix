# latticeBox2D

add_executable(latticeBox2D
  "latticeBox2D.C"
  "writeBasicMesh.C"
  "genericBoundary.C"
  "periodicX.C"
  "periodicY.C"
  "periodicXY.C"  
  )

target_link_libraries(latticeBox2D ${PROJECT_LINK_LIBS} ${MPI_LIBRARIES} ${PROJECT_LINK_LIBS} ${METIS_LIBRARIES})
