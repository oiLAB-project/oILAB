configure_file(bicrystal_2d.txt bicrystal_2d.txt)

# add the executable
add_executable(testGb testGb.cpp)
target_link_libraries(testGb oILAB)

add_test(TestGb testGb)

