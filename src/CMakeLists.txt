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
)
                            
target_link_libraries(oILAB PUBLIC Eigen3::Eigen)

