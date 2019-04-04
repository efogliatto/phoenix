#------------------ Additional calculations ------------------#

add_executable(extraFD "extraFD.C")

target_link_libraries(extraFD ${PROJECT_LINK_LIBS} ${MPI_LIBRARIES} ${PROJECT_LINK_LIBS})
