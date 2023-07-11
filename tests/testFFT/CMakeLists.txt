cmake_minimum_required(VERSION 3.5)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

configure_file(m3D.txt m3D.txt)
configure_file(fft_m3D.txt fft_m3D.txt)
configure_file(m2D.txt m2D.txt)
configure_file(fft_m2D.txt fft_m2D.txt)

# add the executable
add_executable(testFFT testFFT.cpp)
target_link_libraries(testFFT oILAB)

add_test(TestFFT testFFT)

