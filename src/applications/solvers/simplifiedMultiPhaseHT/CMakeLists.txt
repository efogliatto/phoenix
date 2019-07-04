#------------------ Multiphase solver within pseudopotencial framework with heat transfer ------------------#

add_executable(simplifiedMultiPhaseHT "simplifiedMultiPhaseHT.C")

target_link_libraries(simplifiedMultiPhaseHT ${PROJECT_LINK_LIBS} ${MPI_LIBRARIES} ${ARMADILLO_LIBRARIES} ${PROJECT_LINK_LIBS})
