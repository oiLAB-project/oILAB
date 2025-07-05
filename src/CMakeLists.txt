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
                            

target_link_libraries(oILAB PRIVATE Python3::Python)

pybind11_add_module(pyoilab Bindings/module.cpp)
target_link_libraries(pyoilab PRIVATE oILAB)

configure_file(Bindings/example.py example.py)

