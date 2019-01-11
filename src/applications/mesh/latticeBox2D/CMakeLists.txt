#------------------ Lattice Box 2D------------------#


# latticeBox2D

add_executable(latticeBox2D "latticeBox2D.C" # "genericBoundary.c" "periodicX.c" "periodicY.c" "periodicXY.c"
  )

target_link_libraries(latticeBox2D ${PROJECT_LINK_LIBS} ${MPI_LIBRARIES} ${PROJECT_LINK_LIBS})
