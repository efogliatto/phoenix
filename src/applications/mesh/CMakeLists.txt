# latticeBox2D
add_subdirectory (latticeBox2D)

# latticeBox3D
add_subdirectory (latticeBox3D)

# latticeMeshPartition
add_subdirectory (latticeMeshPartitions)

# salomeToPhoenix
add_subdirectory (salomeToPhoenix)



# cartesianMesh

if(NOT CGAL_BASE_DIR STREQUAL "none")

  add_subdirectory (cartesianMesh)

endif()


