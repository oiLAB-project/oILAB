#include <pybind11/pybind11.h>
#include <PyLatticeModule.h>

namespace py = pybind11;

namespace gbLAB {
    PYBIND11_MODULE(gblab_bindings, m) {
        m.doc() = "Python bindings for the oILAB C++ library";
        bind_Lattice<2>(m);
        bind_LatticeVector<2>(m);
        bind_ReciprocalLatticeVector<2>(m);
        bind_LatticeDirection<2>(m);
        bind_Lattice<3>(m);
        bind_LatticeVector<3>(m);
        bind_ReciprocalLatticeVector<3>(m);
        bind_LatticeDirection<3>(m);
    }
}