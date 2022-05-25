include_directories(${gbLABrootDir}/include/Utilities/)
target_sources(${PROJECT_NAME} PRIVATE 
	${CMAKE_CURRENT_SOURCE_DIR}/BiCrystal.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/Lattice.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/LatticeDirection.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/LatticeCore.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/LatticeVector.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/LLL.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/RationalLatticeDirection.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/RationalMatrix.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/ReciprocalLatticeDirection.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/ReciprocalLatticeVector.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/RLLL.cpp
	)

