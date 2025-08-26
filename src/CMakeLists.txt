find_package(pybind11 REQUIRED)
add_library(oILAB SHARED    Lattices/BiCrystal.cpp
                            Lattices/Gb.cpp
                            Lattices/Lattice.cpp 
                            Lattices/LatticeVector.cpp
                            Math/RationalMatrix.cpp
                            Lattices/LatticeCore.cpp
                            Math/RLLL.cpp
                            Lattices/ReciprocalLatticeDirection.cpp
                            Math/LLL.cpp
                            Lattices/LatticeDirection.cpp
                            Lattices/RationalLatticeDirection.cpp
                            Lattices/RationalReciprocalLatticeDirection.cpp
                            Lattices/ReciprocalLatticeVector.cpp
                            Math/Farey.cpp
                            Lattices/GbShifts.cpp
                            Lattices/GbMaterialTensors.cpp)


# Conditionally apply the export property for MSVC on Windows
if(MSVC)
    set_property(TARGET oILAB PROPERTY WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

target_include_directories(oILAB
    PUBLIC
        ${PROJECT_INCLUDE_DIRS}
        ${Python3_INCLUDE_DIRS}
)
                            

target_link_libraries(oILAB 
    PRIVATE Python3::Python 
    PUBLIC Eigen3::Eigen
    PUBLIC OpenMP::OpenMP_CXX)

#pybind11_add_module(pyoilab Bindings/module.cpp)
#target_link_libraries(pyoilab PRIVATE oILAB)

configure_file(Bindings/example.py example.py)
configure_file(Bindings/pyoilab.pyi pyoilab.pyi)

