cmake_minimum_required(VERSION 3.5)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

configure_file(meyer.txt meyer.txt)
configure_file(kurokawa.txt kurokawa.txt)

# add the executable
add_executable(graphene graphene.cpp)
target_link_libraries(graphene oILAB OpenMP::OpenMP_CXX)

