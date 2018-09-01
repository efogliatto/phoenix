#------------------ Unit tests ------------------#


# Test 1. latticeModel
add_executable(test1 "test1.C")
target_link_libraries(test1 ${PROJECT_LINK_LIBS} ${PROJECT_LINK_LIBS})


# Test 2. algebra
add_executable(test2 "test2.C")
target_link_libraries(test2 ${PROJECT_LINK_LIBS} ${PROJECT_LINK_LIBS})
