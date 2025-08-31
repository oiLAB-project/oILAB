# ---------- OpenMP ----------
option(USE_OpenMP "Use OpenMP" ON)
if(USE_OpenMP)
    find_package(OpenMP REQUIRED)
endif()

configure_file(lammps.py lammps.py COPYONLY)
configure_file(Cu_mishin1.eam.alloy Cu_mishin1.eam.alloy COPYONLY)

add_executable(testMonteCarlo testMonteCarlo.cpp)
target_link_libraries(testMonteCarlo
    PRIVATE oILAB
    PRIVATE OpenMP::OpenMP_CXX
)
add_test(TestMonteCarlo testMonteCarlo)

