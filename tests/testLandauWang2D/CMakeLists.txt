configure_file(bicrystal_3d.txt bicrystal_3d.txt)
configure_file(lammps.py lammps.py COPYONLY)
configure_file(Cu_mishin1.eam.alloy Cu_mishin1.eam.alloy COPYONLY)

# add the executable
add_executable(testLandauWang2D testLandauWang2D.cpp)

#target_link_libraries(testLandauWang2D PRIVATE Python3::Python PUBLIC oILAB)
target_link_libraries(testLandauWang2D PUBLIC oILAB)

add_test(TestLandauWang2D testLandauWang2D)

