configure_file(bicrystal_3d.txt bicrystal_3d.txt)
configure_file(lammps.py lammps.py COPYONLY)
configure_file(Cu_mishin1.eam.alloy Cu_mishin1.eam.alloy COPYONLY)

# add the executable
add_executable(testGbMesoState testGbMesoState.cpp)
target_link_libraries(testGbMesoState oILAB)

add_test(TestGbMesoState testGbMesoState)

