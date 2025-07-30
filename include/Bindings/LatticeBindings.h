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
        using PyReciprocalLatticeDirection= PyReciprocalLatticeDirection<dim>;
        using PyLatticeDirection= PyLatticeDirection<dim>;
        using LatticeVector = gbLAB::LatticeVector<dim>;

        using MatrixDimD = Eigen::Matrix<double, dim, dim>;
        using VectorDimD = Eigen::Matrix<double, dim, 1>;

        py::class_<Lattice> cls(m, ("Lattice" + std::to_string(dim) + "D").c_str());
        cls.def(py::init<const MatrixDimD&, const MatrixDimD&>(),
                py::arg("A"), py::arg("Q")=MatrixDimD::Identity())
            .def_readonly("latticeBasis", &Lattice::latticeBasis)
            .def_readonly("reciprocalBasis", &Lattice::reciprocalBasis)
            .def_readonly("F", &Lattice::F)
            .def("interPlanarSpacing", &Lattice::interPlanarSpacing)
            .def("latticeVector", [](const Lattice &lattice, const VectorDimD &p) {
                return PyLatticeVector(lattice.latticeVector(p));
            })
            .def("box",[](const Lattice& lattice, const std::vector<PyLatticeVector>& boxPyLatticeVectors, const std::string& filename){
                std::vector<LatticeVector> boxLatticeVectors;
                for(const auto& v : boxPyLatticeVectors)
                    boxLatticeVectors.push_back(v.lv);
                auto latticeVectors= lattice.box(boxLatticeVectors,filename);

                std::vector<PyLatticeVector> pyLatticeVectors;
                for(const auto& v : latticeVectors)
                    pyLatticeVectors.push_back(PyLatticeVector(v));
                return pyLatticeVectors;
            }, py::arg("boxVectors"),py::arg("filename")="");
        if constexpr(dim==3) {
            cls.def("generateCoincidentLattices",
                 [](const Lattice &lattice, const PyReciprocalLatticeDirection& rd, const double& maxDen, const int& N) {
                     return lattice.generateCoincidentLattices(rd.rld, maxDen, N);
                 }, py::arg("rd"), py::arg("maxDen") = 100, py::arg("N") = 100);
        }
        if constexpr(dim==2) {
            cls.def("generateCoincidentLattices",
                 [](const Lattice& lattice, const double& maxStrain, const double& maxDen, const int& N) {
                     return lattice.generateCoincidentLattices(maxStrain, maxDen, N);
                 }, py::arg("maxStrain"),py::arg("maxDen") = 50, py::arg("N") = 30);
            cls.def("generateCoincidentLattices",
                    [](const Lattice& lattice, const Lattice& underformedLattice, const double& maxStrain, const double& maxDen, const int& N) {
                        return lattice.generateCoincidentLattices(underformedLattice, maxStrain, maxDen, N);
                    }, py::arg("undeformedLattice"), py::arg("maxStrain"),py::arg("maxDen") = 50, py::arg("N") = 30);
        }
        cls.def("latticeDirection",[](const Lattice& self, const VectorDimD& d, const double& tol){
            return PyLatticeDirection(self.latticeDirection(d,tol));
        }, py::arg("cartesianCoordinate"), py::arg("tol")=FLT_EPSILON);
        cls.def("reciprocalLatticeDirection",[](const Lattice& self, const VectorDimD& d, const double& tol){
            return PyReciprocalLatticeDirection(self.reciprocalLatticeDirection(d,tol));
        }, py::arg("cartesianCoordinate"), py::arg("tol")=FLT_EPSILON);
        cls.def("planeParallelLatticeBasis",[](const Lattice& self, const PyReciprocalLatticeDirection& l, const bool& useRLLL){
            auto latticeBasis2D= self.planeParallelLatticeBasis(l.rld, useRLLL);
            std::vector<PyLatticeDirection> pyLatticeBasis2D;
            for(const auto elem : latticeBasis2D)
                pyLatticeBasis2D.push_back(elem);
            return pyLatticeBasis2D;
        }, py::arg("reciprocalLatticeDirection"), py::arg("useRLLL")=false);
    }
}
#endif //OILAB_LATTICE_BINDINGS_H
