cmake_minimum_required(VERSION 3.20)
project(oILAB VERSION 1.0 LANGUAGES C CXX)

# ---------- C++ standard ----------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)

if(MSVC)
    add_compile_options(/bigobj /MP)
endif()

# ---------- Dependencies ----------
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# ---------- Subdirectories ----------
add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(userExamples)

# ---------- Testing (optional) ----------
option(ENABLE_TESTING "Build tests" OFF)
if(ENABLE_TESTING)
    add_subdirectory(tests)
endif()

 # ---------- Python bindings (optional) ----------
option(BUILD_PYBINDINGS "Build Python bindings (pyOILAB)" ON)
if(BUILD_PYBINDINGS)
    add_subdirectory(bindings)
endif()

# ---------- Installation ----------
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Install header files (preserves subdirectory structure under oILAB/)
install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/oILAB
    FILES_MATCHING PATTERN "*.h"
)

# Generate package config files
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/oILABConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/oILABConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/oILAB
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/oILABConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

# Install CMake config files
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/oILABConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/oILABConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/oILAB
)

# Export targets
install(EXPORT oILABTargets
    FILE oILABTargets.cmake
    NAMESPACE oILAB::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/oILAB
)
