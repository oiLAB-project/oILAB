cmake_minimum_required(VERSION 3.5)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# add the executable
add_executable(testBezout testBezout.cpp)
target_link_libraries(testBezout oILAB)

add_test(TestBezout testBezout)

