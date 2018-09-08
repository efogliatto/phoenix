#------------------ Unit tests ------------------#


# Test 1. latticeModel
add_executable(test1 "test1.C")
target_link_libraries(test1 ${PROJECT_LINK_LIBS} ${PROJECT_LINK_LIBS})


# Test 2. algebra
add_executable(test2 "test2.C")
target_link_libraries(test2 ${PROJECT_LINK_LIBS} ${PROJECT_LINK_LIBS})


# Test 3. MAtrix - vector multiplication efficiency in c
add_executable(test3 "test3.c")
target_link_libraries(test3 ${PROJECT_LINK_LIBS} ${PROJECT_LINK_LIBS})


# Test 4. Dictionary
add_executable(test4 "test4.C")
target_link_libraries(test4 ${PROJECT_LINK_LIBS} ${PROJECT_LINK_LIBS})


# Test 5. mpiInfo
add_executable(test5 "test5.C")
target_link_libraries(test5 ${PROJECT_LINK_LIBS} ${MPI_LIBRARIES} ${PROJECT_LINK_LIBS})
