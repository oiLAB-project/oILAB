cmake_minimum_required(VERSION 3.5)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# add the executable
add_executable(testLattice testLattice.cpp)
target_link_libraries(testLattice oILAB)

add_test(TestLattice testLattice)

