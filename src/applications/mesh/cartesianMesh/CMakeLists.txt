# cartesianMesh

add_executable(cartesianMesh
  "cartesianMesh.C"
  "STLToPolyhedron.C"
  )

target_link_libraries(cartesianMesh ${PROJECT_LINK_LIBS} ${MPI_LIBRARIES} ${PROJECT_LINK_LIBS} CGAL::CGAL)
