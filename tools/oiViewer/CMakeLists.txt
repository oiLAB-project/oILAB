cmake_minimum_required(VERSION 3.1.0)
project(oiViewer VERSION 1.0.0 LANGUAGES CXX)

set(EigenRootDir "/usr/local/include")
set(CMAKE_PREFIX_PATH "/Applications/Qt/6.1.3/macos/lib/cmake/")
set(gbLABrootDir "${CMAKE_CURRENT_SOURCE_DIR}/../../")
set(CMAKE_CXX_FLAGS "-Wall -Wextra -O3")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

add_executable(${PROJECT_NAME} oiViewer.cxx)

find_package(Qt6 COMPONENTS Widgets OpenGL REQUIRED)
message (STATUS "Qt6_VERSION: ${Qt6_VERSION}")
if (NOT Qt6_FOUND)
  message("Skipping DDvtk: ${Qt6_NOT_FOUND_MESSAGE}")
  return()
endif()

find_package(VTK REQUIRED)
message (STATUS "VTK_VERSION: ${VTK_VERSION}")
if (NOT VTK_FOUND)
  message("Skipping DDvtk: ${VTK_NOT_FOUND_MESSAGE}")
  return()
endif()

include_directories(${EigenRootDir})
include_directories(${gbLABrootDir}/include/IO/)
include_directories(${gbLABrootDir}/include/Utilities/)
include_directories(${gbLABrootDir}/include/Math/)
include_directories(${gbLABrootDir}/include/Lattices/)
include_directories(${gbLABrootDir}/include/Visualization/)

add_subdirectory(${gbLABrootDir}/src/Lattices Lattices)
add_subdirectory(${gbLABrootDir}/src/Visualization Visualization)


target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES} Qt6::Widgets Qt6::OpenGL Qt6::Core)
vtk_module_autoinit(TARGETS oiViewer MODULES ${VTK_LIBRARIES} )
