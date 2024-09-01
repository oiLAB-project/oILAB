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
                            GbContinuum.cpp
                            GbMaterialTensors.cpp
                            GbMesoState.cpp
                            GbMesoStateEnsemble.cpp
                            ../MonteCarlo/MonteCarlo.cpp
                            ../MonteCarlo/EvolutionAlgorithm.cpp
                            ../MonteCarlo/CanonicalTP.cpp
        ../MonteCarlo/LandauWangTP.cpp)

target_link_libraries(oILAB PRIVATE Python3::Python)


