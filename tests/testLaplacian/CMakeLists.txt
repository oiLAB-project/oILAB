cmake_minimum_required(VERSION 3.5)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# add the executable
add_executable(testLaplacian testLaplacian.cpp)
target_link_libraries(testLaplacian oILAB)

add_test(TestLaplacian testLaplacian)

