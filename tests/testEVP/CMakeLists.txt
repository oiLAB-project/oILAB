cmake_minimum_required(VERSION 3.5)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# add the executable
add_executable(testEVP testEVP.cpp)
target_link_libraries(testEVP oILAB)

add_test(TestEVP testEVP)

