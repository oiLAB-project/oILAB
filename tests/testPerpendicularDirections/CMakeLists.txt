cmake_minimum_required(VERSION 3.5)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# add the executable
add_executable(testPerpendicularDirections testPerpendicularDirections.cpp)
target_link_libraries(testPerpendicularDirections oILAB)

add_test(TestPerpendicularDirections testPerpendicularDirections)

