include_directories(${gbLABrootDir}/include/Utilities/)
include_directories(${gbLABrootDir}/include/Lattices/)
include_directories(${gbLABrootDir}/include/IO/)
include_directories(${gbLABrootDir}/include/Math/)
include_directories(${gbLABrootDir}/include/Visualization/)


target_sources(${PROJECT_NAME} PRIVATE 
	${gbLABrootDir}/include/Visualization/oiViewerMainWindow.h
	${CMAKE_CURRENT_SOURCE_DIR}/oiViewerMainWindow.cpp
	)
