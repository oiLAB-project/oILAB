#include <pybind11/pybind11.h>

namespace py = pybind11;

void bind_Lattice(py::module_ &);

PYBIND11_MODULE(gblab_bindings, m) {
    m.doc() = "Python bindings for the oILAB C++ library";
    bind_Lattice(m);
}

