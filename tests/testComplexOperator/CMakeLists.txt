cmake_minimum_required(VERSION 3.5)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

configure_file(evalues.txt evalues.txt)
configure_file(complexMatrix.txt complexMatrix.txt)
configure_file(realMatrix.txt realMatrix.txt)

# add the executable
add_executable(testComplexOperator testComplexOperator.cpp)
target_link_libraries(testComplexOperator oILAB)

add_test(TestComplexOperator testComplexOperator)

