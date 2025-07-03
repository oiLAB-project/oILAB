//
// Created by Nikhil Chandra Admal on 7/3/25.
//

#ifndef OILAB_LATTICE_BINDINGS_H
#define OILAB_LATTICE_BINDINGS_H

#include <pybind11/pybind11.h>
#include <LatticeModule.h>
#include <LatticeVectorBindings.h>
namespace py = pybind11;

namespace gbLAB {
    template<int dim>
    void bind_Lattice(py::module_ &m) {
        using Lattice = gbLAB::Lattice<dim>;
        using PyLatticeVector = PyLatticeVector<dim>;
        using MatrixDimD = Eigen::Matrix<double, dim, dim>;
        using VectorDimD = Eigen::Matrix<double, dim, 1>;

        py::class_<Lattice>(m, ("Lattice" + std::to_string(dim) + "D").c_str())
                .def(py::init<const MatrixDimD&, const MatrixDimD&>(),
                        py::arg("A"), py::arg("Q")=MatrixDimD::Identity())
                .def_readonly("latticeBasis", &Lattice::latticeBasis)
                .def_readonly("reciprocalBasis", &Lattice::reciprocalBasis)
                .def_readonly("F", &Lattice::F)
                .def("interPlanarSpacing", &Lattice::interPlanarSpacing)
                .def("latticeVector", [](const Lattice &lattice, const VectorDimD &p) {
                    return PyLatticeVector(lattice.latticeVector(p));
                });
    }
}
#endif //OILAB_LATTICE_BINDINGS_H
