configure_file(bicrystal_2d.txt bicrystal_2d.txt)

# add the executable
add_executable(testMoire testMoire.cpp)
target_link_libraries(testMoire oILAB)

add_test(TestMoire testMoire)

