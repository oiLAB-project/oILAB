# add the executable
add_executable(testCCUM testCCUM.cpp)
target_link_libraries(testCCUM oILAB)

add_test(TestCCUM testCCUM)

