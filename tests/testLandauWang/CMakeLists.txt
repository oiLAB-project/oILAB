cmake_minimum_required(VERSION 3.5)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

configure_file(bicrystal_3d.txt bicrystal_3d.txt)
configure_file(lammps.py lammps.py COPYONLY)
configure_file(Cu_mishin1.eam.alloy Cu_mishin1.eam.alloy COPYONLY)

# add the executable
add_executable(testLandauWang testLandauWang.cpp)

#target_link_libraries(testLandauWang PRIVATE Python3::Python PUBLIC oILAB)
target_link_libraries(testLandauWang PUBLIC oILAB)

add_test(TestLandauWang testLandauWang)

