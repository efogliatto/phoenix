#------------------ Multiphase solver within pseudopotencial framework with heat transfer ------------------#

add_executable(multiPhasePPHT "multiPhasePPHT.C")

target_link_libraries(multiPhasePPHT ${PROJECT_LINK_LIBS} ${MPI_LIBRARIES} ${PROJECT_LINK_LIBS})
