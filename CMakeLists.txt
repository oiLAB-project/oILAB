cmake_minimum_required(VERSION 3.20)
project(oILAB VERSION 1.0 LANGUAGES C CXX)

# ---------- C++ standard ----------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Debug)

if(MSVC)
    add_compile_options(/bigobj /MP)
endif()

# ---------- Dependencies ----------
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# ---------- Project include directories ----------
set(PROJECT_INCLUDE_DIRS
    ${PROJECT_SOURCE_DIR}/include/IO
    ${PROJECT_SOURCE_DIR}/include/Lattices
    ${PROJECT_SOURCE_DIR}/include/MonteCarlo
    ${PROJECT_SOURCE_DIR}/include/Math
    ${PROJECT_SOURCE_DIR}/include/Utilities
    ${PROJECT_SOURCE_DIR}/include/Visualization
    ${PROJECT_SOURCE_DIR}/include/Bindings
)

# ---------- Subdirectories ----------
add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(examples)
add_subdirectory(userExamples)

 # ---------- Python bindings (optional) ----------
option(BUILD_PYBINDINGS "Build Python bindings (pyOILAB)" ON)
if(BUILD_PYBINDINGS)
    add_subdirectory(bindings)
endif()
