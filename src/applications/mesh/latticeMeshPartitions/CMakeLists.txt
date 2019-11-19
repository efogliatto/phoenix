# latticeMeshPartition

add_executable(latticeMeshPartitions "latticeMeshPartitions.C" # "standardDecomp.C" "kmetisDecomp.C" "localIndexing.C"
  )

target_link_libraries(latticeMeshPartitions ${PROJECT_LINK_LIBS} ${MPI_LIBRARIES} ${PROJECT_LINK_LIBS})
