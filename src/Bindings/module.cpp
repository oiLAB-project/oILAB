#include <pybind11/pybind11.h>
#include <PyLatticeModule.h>

namespace py = pybind11;

namespace pyoilab{
    PYBIND11_MODULE(pyoilab, m) {
        m.doc() = "Python bindings for the oILAB C++ library";
        // Dimension 2
        bind_Lattice<2>(m);
        bind_LatticeVector<2>(m);
        bind_LatticeDirection<2>(m);
        bind_ReciprocalLatticeVector<2>(m);
        bind_ReciprocalLatticeDirection<2>(m);

        // Dimension 3
        bind_Lattice<3>(m);
        bind_LatticeVector<3>(m);
        bind_LatticeDirection<3>(m);
        bind_ReciprocalLatticeVector<3>(m);
        bind_ReciprocalLatticeDirection<3>(m);
    }
}