cmake_minimum_required(VERSION 3.0)
project(bicrystal CXX)

set(EigenRootDir "/usr/local/include")
set(gbLABrootDir "${CMAKE_CURRENT_SOURCE_DIR}/../..")
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
message("CMAKE_SYSTEM_NAME IS: ${CMAKE_SYSTEM_NAME}")
if(${CMAKE_SYSTEM_NAME} MATCHES Darwin)
#set(CMAKE_CXX_COMPILER "/opt/local/bin/g++")
	set(CMAKE_CXX_FLAGS "-Wall -Wextra -O3 -Ofast")
elseif(${CMAKE_SYSTEM_NAME} MATCHES Linux)
	set(CMAKE_CXX_FLAGS "-Wall -Wextra -march=native")
endif()

add_executable(bicrystal main.cpp)

include_directories(/usr/local/include)
include_directories(${gbLABrootDir}/include/IO/)
include_directories(${gbLABrootDir}/include/Utilities/)
include_directories(${gbLABrootDir}/include/Math/)
include_directories(${gbLABrootDir}/include/Lattices/)

add_subdirectory(${gbLABrootDir}/src/Lattices Lattices)





