# latticeMeshPartition

add_executable(salomeToPhoenix
  "salomeToPhoenix.C"
  "readBasicMesh.C"
  )

target_link_libraries(salomeToPhoenix ${PROJECT_LINK_LIBS} ${MPI_LIBRARIES} ${PROJECT_LINK_LIBS})
