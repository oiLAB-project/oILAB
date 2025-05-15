cmake_minimum_required(VERSION 3.5)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

configure_file(cubic_lattice.txt cubic_lattice.txt)

# add the executable
add_executable(testSNF_cubic testSNF_cubic.cpp)
target_link_libraries(testSNF_cubic oILAB)

add_test(TestSNF_cubic testSNF_cubic)

