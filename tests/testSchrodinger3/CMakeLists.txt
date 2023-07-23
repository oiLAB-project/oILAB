cmake_minimum_required(VERSION 3.5)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

configure_file(meyer.txt meyer.txt)
configure_file(kurokawa.txt kurokawa.txt)

# add the executable
add_executable(testSchrodinger3 testSchrodinger3.cpp)
target_link_libraries(testSchrodinger3 oILAB OpenMP::OpenMP_CXX)
add_test(TestSchrodinger3 testSchrodinger3)

