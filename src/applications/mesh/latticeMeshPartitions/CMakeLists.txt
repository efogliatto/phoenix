# latticeMeshPartition

add_executable(latticeMeshPartitions
  "latticeMeshPartitions.C"
  "readBasicMesh.C"
  "standardDecomp.C"
  "kmetisDecomp.C"
  "mpmetisDecomp.C"
  "mlevelKmetisDecomp.C"
  "localIndexing.C"
  "writeLatticeMesh.C"
  "writeMeshToEnsight.C"
  "updateCaseFile.C"
  "computeVirtualNodes.C"
  )

target_link_libraries(latticeMeshPartitions ${PROJECT_LINK_LIBS} ${MPI_LIBRARIES} ${PROJECT_LINK_LIBS})
