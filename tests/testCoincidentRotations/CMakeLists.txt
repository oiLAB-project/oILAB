configure_file(bicrystal_3d.txt bicrystal_3d.txt)

# add the executable
add_executable(testCoincidentRotations testCoincidentRotations.cpp)
target_link_libraries(testCoincidentRotations oILAB)

add_test(TestCoincidentRotations testCoincidentRotations)

