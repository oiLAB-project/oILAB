//
// Created by Nikhil Chandra Admal on 7/3/25.
//

#ifndef OILAB_LATTICE_BINDINGS_H
#define OILAB_LATTICE_BINDINGS_H

#include <pybind11/pybind11.h>
#include <LatticeModule.h>
#include <pybind11/stl.h>
namespace py = pybind11;

namespace pyoilab {
    template<int dim>
    void bind_Lattice(py::module_ &m) {
        using Lattice = gbLAB::Lattice<dim>;
        using PyLatticeVector = PyLatticeVector<dim>;
        using LatticeVector = gbLAB::LatticeVector<dim>;

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
                })
                .def("box",[](const Lattice& lattice, const std::vector<PyLatticeVector>& boxPyLatticeVectors, const std::string& filename=""){
                    std::vector<LatticeVector> boxLatticeVectors;
                    for(const auto& v : boxPyLatticeVectors)
                        boxLatticeVectors.push_back(v.lv);
                    auto latticeVectors= lattice.box(boxLatticeVectors,filename);

                    std::vector<PyLatticeVector> pyLatticeVectors;
                    for(const auto& v : latticeVectors)
                        pyLatticeVectors.push_back(PyLatticeVector(v));
                    return pyLatticeVectors;
                }, py::arg("boxVectors"),py::arg("filename")="");

                /*
        template<int dm=dim>
        typename std::enable_if<dm==3,std::vector<LatticeVector<dim>>>::type
        box(const std::vector<LatticeVector<dim>>& boxVectors, const std::string& filename= "") const;
                 */
    }
}
#endif //OILAB_LATTICE_BINDINGS_H
