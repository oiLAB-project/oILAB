# add the executable
add_executable(testBezout testBezout.cpp)
target_link_libraries(testBezout oILAB)

add_test(TestBezout testBezout)

