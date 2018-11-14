#------------------ Temperature gradient on non periodic walls ------------------#

add_executable(wallHeatFlux "wallHeatFlux.C")

target_link_libraries(wallHeatFlux ${PROJECT_LINK_LIBS} ${MPI_LIBRARIES} ${PROJECT_LINK_LIBS})
