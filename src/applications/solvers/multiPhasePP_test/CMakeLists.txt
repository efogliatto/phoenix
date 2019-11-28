#------------------ Multiphase solver within pseudopotencial framework ------------------#

add_executable(multiPhasePP_test "multiPhasePP_test.C")

target_link_libraries(multiPhasePP_test ${PROJECT_LINK_LIBS} ${MPI_LIBRARIES} ${PROJECT_LINK_LIBS})
