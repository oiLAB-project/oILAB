cmake_minimum_required(VERSION 3.5)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

configure_file(bicrystal_3d.txt bicrystal_3d.txt)

# add the executable
add_executable(testMonteCarlo testMonteCarlo.cpp)
target_link_libraries(testMonteCarlo oILAB)

add_test(TestMonteCarlo testMonteCarlo)

