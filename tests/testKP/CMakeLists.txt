cmake_minimum_required(VERSION 3.5)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# add the executable
add_executable(testKP testKP.cpp)
target_link_libraries(testKP oILAB OpenMP::OpenMP_CXX)
add_test(TestKP testKP)

