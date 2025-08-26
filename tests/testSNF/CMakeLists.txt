configure_file(bicrystal_2d.txt bicrystal_2d.txt)

# add the executable
add_executable(testSNF testSNF.cpp)
target_link_libraries(testSNF oILAB)

add_test(TestSNF testSNF)

