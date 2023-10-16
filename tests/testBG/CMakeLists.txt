cmake_minimum_required(VERSION 3.5)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

configure_file(meyer.txt meyer.txt)

# add the executable
add_executable(testBG testBG.cpp)
target_link_libraries(testBG oILAB OpenMP::OpenMP_CXX)
add_test(TestBG testBG)

