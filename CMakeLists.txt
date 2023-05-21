cmake_minimum_required(VERSION 3.5)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set the project name
project(oILAB VERSION 1.0)

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
