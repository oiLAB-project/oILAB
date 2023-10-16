cmake_minimum_required(VERSION 3.5)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# set the project name
project(oILAB VERSION 1.0)
set(CMAKE_CXX_COMPILER g++-12)
set(CMAKE_C_COMPILER gcc-12)
#set(CMAKE_CXX_COMPILER g++-mp-12)
#set(CMAKE_C_COMPILER gcc-mp-12)

OPTION (USE_OpenMP "Use OpenMP to enable <omp.h>" ON)
if (USE_OpenMP)
    set(OpenMP_C_FLAGS "-fopenmp")
    set(OpenMP_CXX_FLAGS "-fopenmp")
    set(OpenMP_C_LIB_NAMES "omp")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY omp)
    find_package(OpenMP REQUIRED)
endif(USE_OpenMP)

if(OPENMP_FOUND)
    set(OPENMP_INCLUDES "/usr/lib/llvm-14/include")
    set(OPENMP_LIBRARIES "/usr/lib/llvm-14/lib")
#    set(OPENMP_INCLUDES "/opt/Homebrew/Cellar/libomp/16.0.6/include")
#    set(OPENMP_LIBRARIES "/opt/Homebrew/Cellar/libomp/16.0.6/lib")
    include_directories("${OPENMP_INCLUDES}")
    link_directories("${OPENMP_LIBRARIES}")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OPENMP_FOUND)


# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include/IO)
include_directories(${PROJECT_SOURCE_DIR}/include/Lattices)
include_directories(${PROJECT_SOURCE_DIR}/include/Math)
include_directories(${PROJECT_SOURCE_DIR}/include/Utilities)
include_directories(${PROJECT_SOURCE_DIR}/include/Visualization)
include_directories(${PROJECT_SOURCE_DIR}/include/eigen-3.4.0)
include_directories(${PROJECT_SOURCE_DIR}/include/eigen-3.4.0)
add_subdirectory(src/Lattices)
add_subdirectory(tests)
add_subdirectory(examples)
