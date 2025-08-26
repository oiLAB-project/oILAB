configure_file(bicrystal_3d.txt bicrystal_3d.txt)

# add the executable
add_executable(testMesoStates testMesoStates.cpp)
target_link_libraries(testMesoStates oILAB)

add_test(TestMesoStates testMesoStates)

