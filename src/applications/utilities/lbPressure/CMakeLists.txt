#------------------ Pressure calculation for pseudopotential multiphase models ------------------#

add_executable(lbPressure "lbPressure.C")

target_link_libraries(lbPressure ${PROJECT_LINK_LIBS} ${MPI_LIBRARIES} ${PROJECT_LINK_LIBS})
