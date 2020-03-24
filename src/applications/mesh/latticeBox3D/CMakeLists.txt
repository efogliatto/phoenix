# latticeBox3D

add_executable(latticeBox3D
  "latticeBox3D.C"
  "writeBasicMesh.C"
  "genericBoundary3D.C"
  "periodicXY3D.C"
  "periodicXYZ3D.C"
  "periodicZ3D.C"  
  )

target_link_libraries(latticeBox3D ${PROJECT_LINK_LIBS} ${MPI_LIBRARIES} ${PROJECT_LINK_LIBS} ${METIS_LIBRARIES})
