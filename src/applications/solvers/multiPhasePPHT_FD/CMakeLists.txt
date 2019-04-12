#------------------ Multiphase solver within pseudopotencial framework with heat transfer ------------------#

add_executable(multiPhasePPHT_FD "multiPhasePPHT_FD.C")

target_link_libraries(multiPhasePPHT_FD ${PROJECT_LINK_LIBS} ${MPI_LIBRARIES} ${PROJECT_LINK_LIBS})
