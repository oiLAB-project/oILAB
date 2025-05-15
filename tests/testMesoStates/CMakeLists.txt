cmake_minimum_required(VERSION 3.5)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

configure_file(bicrystal_3d.txt bicrystal_3d.txt)

# add the executable
add_executable(testMesoStates testMesoStates.cpp)
target_link_libraries(testMesoStates oILAB)

add_test(TestMesoStates testMesoStates)

