cmake_minimum_required(VERSION 3.5)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

configure_file(lammps.py lammps.py COPYONLY)
configure_file(Cu_mishin1.eam.alloy Cu_mishin1.eam.alloy COPYONLY)

# add the executable
add_executable(testMonteCarlo testMonteCarlo.cpp)

target_link_libraries(testMonteCarlo PRIVATE Python3::Python PUBLIC oILAB)
#target_link_libraries(testMonteCarlo PUBLIC oILAB)

add_test(TestMonteCarlo testMonteCarlo)

