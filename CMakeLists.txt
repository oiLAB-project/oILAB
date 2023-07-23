cmake_minimum_required(VERSION 3.5)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# set the project name
project(oILAB VERSION 1.0)

OPTION (USE_OpenMP "Use OpenMP to enamble <omp.h>" ON)
if (APPLE AND USE_OpenMP)
    if(CMAKE_C_COMPILER_ID MATCHES "Clang\$")
        set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp")
        set(OpenMP_C_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY omp)
    endif()

    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang\$")
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
        set(OpenMP_CXX_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY omp)
    endif()
endif()

if (USE_OpenMP)
    find_package(OpenMP REQUIRED)
endif(USE_OpenMP)

if(OPENMP_FOUND)
    set(OPENMP_INCLUDES "/opt/Homebrew/Cellar/libomp/16.0.6/include")
    set(OPENMP_LIBRARIES "/opt/Homebrew/Cellar/libomp/16.0.6/lib")
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
