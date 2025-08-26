# add the executable
add_executable(testLattice testLattice.cpp)
target_link_libraries(testLattice oILAB)

add_test(TestLattice testLattice)

