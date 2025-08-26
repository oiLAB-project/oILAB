configure_file(Grimmer1985.txt Grimmer1985.txt)
configure_file(bicrystal_2d.txt bicrystal_2d.txt)
configure_file(bicrystal_3d.txt bicrystal_3d.txt)

# add the executable
add_executable(bicrystal main.cpp)
target_link_libraries(bicrystal oILAB)

