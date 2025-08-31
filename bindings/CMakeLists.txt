find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
find_package(pybind11 REQUIRED)

pybind11_add_module(pyoilab module.cpp)
target_include_directories(pyoilab
    PUBLIC ${Python3_INCLUDE_DIRS}
)
target_link_libraries(pyoilab
    PRIVATE Python3::Python 
    PRIVATE oILAB)

configure_file(example.py example.py)
configure_file(pyoilab.pyi pyoilab.pyi)

