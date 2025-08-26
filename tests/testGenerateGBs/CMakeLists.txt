configure_file(bicrystal_3d.txt bicrystal_3d.txt)

# add the executable
add_executable(testGenerateGBs testGenerateGBs.cpp)
target_link_libraries(testGenerateGBs oILAB)

add_test(TestGenerateGBs testGenerateGBs)

