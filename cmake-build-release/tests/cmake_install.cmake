# Install script for directory: /Users/Nikhil/Documents/Academic/Software/oILAB/tests

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/objdump")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/Nikhil/Documents/Academic/Software/oILAB/cmake-build-release/tests/testSNF/cmake_install.cmake")
  include("/Users/Nikhil/Documents/Academic/Software/oILAB/cmake-build-release/tests/testSNF_cubic/cmake_install.cmake")
  include("/Users/Nikhil/Documents/Academic/Software/oILAB/cmake-build-release/tests/testPerpendicularDirections/cmake_install.cmake")
  include("/Users/Nikhil/Documents/Academic/Software/oILAB/cmake-build-release/tests/testPlaneParallelLatticeDirections/cmake_install.cmake")
  include("/Users/Nikhil/Documents/Academic/Software/oILAB/cmake-build-release/tests/testBezout/cmake_install.cmake")
  include("/Users/Nikhil/Documents/Academic/Software/oILAB/cmake-build-release/tests/testCCUM/cmake_install.cmake")
  include("/Users/Nikhil/Documents/Academic/Software/oILAB/cmake-build-release/tests/testCoincidentRotations/cmake_install.cmake")
  include("/Users/Nikhil/Documents/Academic/Software/oILAB/cmake-build-release/tests/testGb/cmake_install.cmake")
  include("/Users/Nikhil/Documents/Academic/Software/oILAB/cmake-build-release/tests/testLattice/cmake_install.cmake")
  include("/Users/Nikhil/Documents/Academic/Software/oILAB/cmake-build-release/tests/testGenerateGBs/cmake_install.cmake")
  include("/Users/Nikhil/Documents/Academic/Software/oILAB/cmake-build-release/tests/testMoire/cmake_install.cmake")
  include("/Users/Nikhil/Documents/Academic/Software/oILAB/cmake-build-release/tests/testSpectra/cmake_install.cmake")
  include("/Users/Nikhil/Documents/Academic/Software/oILAB/cmake-build-release/tests/testFFT/cmake_install.cmake")
  include("/Users/Nikhil/Documents/Academic/Software/oILAB/cmake-build-release/tests/testEVP/cmake_install.cmake")
  include("/Users/Nikhil/Documents/Academic/Software/oILAB/cmake-build-release/tests/testDiff/cmake_install.cmake")
  include("/Users/Nikhil/Documents/Academic/Software/oILAB/cmake-build-release/tests/testPeriodicFunction/cmake_install.cmake")
  include("/Users/Nikhil/Documents/Academic/Software/oILAB/cmake-build-release/tests/testSchrodinger3/cmake_install.cmake")
  include("/Users/Nikhil/Documents/Academic/Software/oILAB/cmake-build-release/tests/testComplexOperator/cmake_install.cmake")

endif()

