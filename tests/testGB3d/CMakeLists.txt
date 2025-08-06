cmake_minimum_required(VERSION 3.5)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

configure_file(bicrystal_2d.txt bicrystal_2d.txt)

# add the executable
add_executable(testGb3d testGb3d.cpp)
target_link_libraries(testGb3d oILAB)

add_test(TestGb3d testGb3d)

