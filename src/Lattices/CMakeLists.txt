add_library(oILAB SHARED    BiCrystal.cpp
                            Gb.cpp
                            Lattice.cpp 
                            LatticeVector.cpp
                            ../Math/RationalMatrix.cpp
                            LatticeCore.cpp
                            ../Math/RLLL.cpp
                            ReciprocalLatticeDirection.cpp
                            ../Math/LLL.cpp
                            LatticeDirection.cpp
                            RationalLatticeDirection.cpp
                            RationalReciprocalLatticeDirection.cpp
                            ReciprocalLatticeVector.cpp
                            ../Math/Farey.cpp
                            GbShifts.cpp
                            GbMaterialTensors.cpp)

target_link_libraries(oILAB PRIVATE Python3::Python)


