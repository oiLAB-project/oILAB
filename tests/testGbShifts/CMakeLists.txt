configure_file(bicrystal_3d.txt bicrystal_3d.txt)

# add the executable
add_executable(testGbShifts testGbShifts.cpp)
target_link_libraries(testGbShifts oILAB)

add_test(TestGbShifts testGbShifts)

