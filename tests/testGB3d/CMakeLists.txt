configure_file(bicrystal_2d.txt bicrystal_2d.txt)

# add the executable
add_executable(testGb3d testGb3d.cpp)
target_link_libraries(testGb3d oILAB)

add_test(TestGb3d testGb3d)

