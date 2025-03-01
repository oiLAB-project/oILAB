cmake_minimum_required(VERSION 3.5)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

enable_testing()
add_subdirectory(testSNF)
add_subdirectory(testSNF_cubic)
add_subdirectory(testPerpendicularDirections)
add_subdirectory(testPlaneParallelLatticeDirections)
add_subdirectory(testBezout)
add_subdirectory(testCCUM)
add_subdirectory(testCoincidentRotations)
add_subdirectory(testGb)
add_subdirectory(testGB3d)
add_subdirectory(testLattice)
add_subdirectory(testGenerateGBs)
add_subdirectory(testMoire)
