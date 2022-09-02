cmake_minimum_required(VERSION 3.5)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

configure_file(bicrystal_2d.txt bicrystal_2d.txt)

# add the executable
add_executable(testSNF testSNF.cpp)
target_link_libraries(testSNF oILAB)

add_test(TestSNF testSNF)

