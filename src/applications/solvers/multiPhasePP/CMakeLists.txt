#------------------ Multiphase solver within pseudopotencial framework ------------------#

add_executable(multiPhasePP "multiPhasePP.C")

target_link_libraries(multiPhasePP ${PROJECT_LINK_LIBS} ${MPI_LIBRARIES} ${PROJECT_LINK_LIBS})
