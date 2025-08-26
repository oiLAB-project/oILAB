cmake_minimum_required(VERSION 3.20)
project(oILAB VERSION 1.0 LANGUAGES C CXX)

# ---------- Compiler ----------
# Can be overridden from terminal via -D or environment variables
# Example: -DCMAKE_C_COMPILER=/opt/homebrew/opt/llvm/bin/clang
#          -DCMAKE_CXX_COMPILER=/opt/homebrew/opt/llvm/bin/clang++
set(CMAKE_C_COMPILER "" CACHE STRING "C compiler")
set(CMAKE_CXX_COMPILER "" CACHE STRING "C++ compiler")

# ---------- C++ standard ----------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Debug)

# ---------- OpenMP ----------
# Can be overridden from terminal via -D or environment variables
# Example:  -DUSE_OpenMP=ON
#           -DOpenMP_C_FLAGS="-fopenmp"
#           -DOpenMP_CXX_FLAGS="-fopenmp"
#           -DOpenMP_C_LIB_NAMES="omp"
#           -DOpenMP_CXX_LIB_NAMES="omp"
#           -DOpenMP_omp_LIBRARY="/opt/homebrew/opt/llvm/lib/libomp.dylib"
option(USE_OpenMP "Use OpenMP" ON)
if(USE_OpenMP)
    find_package(OpenMP REQUIRED)
endif()

if(MSVC)
    add_compile_options(/bigobj /MP)
endif()

# ---------- Dependencies ----------
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)

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

