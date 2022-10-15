cmake_minimum_required(VERSION 3.5)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

configure_file(bicrystal_3d.txt bicrystal_3d.txt)

# add the executable
add_executable(testCoincidentRotations testCoincidentRotations.cpp)
target_link_libraries(testCoincidentRotations oILAB)

add_test(TestCoincidentRotations testCoincidentRotations)

