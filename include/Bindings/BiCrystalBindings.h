//
// Created by Nikhil Chandra Admal on 7/9/25.
//

#ifndef OILAB_BICRYSTAL_BINDINGS_H
#define OILAB_BICRYSTAL_BINDINGS_H

#include <pybind11/pybind11.h>
#include <LatticeModule.h>
#include <pybind11/stl.h>
namespace py = pybind11;

namespace pyoilab {
    template<int dim>
    void bind_BiCrystal(py::module_ &m) {
        using PyLatticeVector = PyLatticeVector<dim>;
        using BiCrystal = gbLAB::BiCrystal<dim>;
        using Lattice = gbLAB::Lattice<dim>;
        using LatticeVector = gbLAB::LatticeVector<dim>;

        using MatrixDimD = Eigen::Matrix<double, dim, dim>;
        using VectorDimD = Eigen::Matrix<double, dim, 1>;

        py::class_<BiCrystal> cls(m, ("BiCrystal" + std::to_string(dim) + "D").c_str());
        cls.def(py::init<const Lattice&, const Lattice&, const bool&>(),
                py::arg("A"), py::arg("B"), py::arg("useRLLL")=false);
        cls.def("A",[](const BiCrystal& self) -> const Lattice& {
            return self.A;
        }, py::return_value_policy::reference_internal);
        cls.def("B",[](const BiCrystal& self) -> const Lattice& {
            return self.B;
        }, py::return_value_policy::reference_internal);
        cls.def_readonly("sigma",&BiCrystal::sigma);
        cls.def_readonly("csl",&BiCrystal::csl);
        cls.def_readonly("dscl",&BiCrystal::dscl);
        cls.def("box", [](const BiCrystal& self,
                          std::vector<PyLatticeVector>& boxPyLatticeVectors,
                          const double& orthogonality,
                          const int& dsclFactor,
                          std::string filename,
                          bool orient){
            std::vector<LatticeVector> boxLatticeVectors;
            for(const auto& v : boxPyLatticeVectors)
                boxLatticeVectors.push_back(v.lv);
            auto latticeVectors= self.box(boxLatticeVectors,
                                          orthogonality,
                                          dsclFactor,
                                          filename,
                                          orient);

            std::vector<PyLatticeVector> pyLatticeVectors;
            for(const auto& v : latticeVectors)
                pyLatticeVectors.push_back(PyLatticeVector(v));
            return pyLatticeVectors;
        }, py::arg("boxVectors"), py::arg("orthogonality"), py::arg("dsclFactor"), py::arg("filename")="", py::arg("orient")=false);
        cls.def("getLatticeDirectionInC",[](const BiCrystal& self, const PyLatticeVector& v){
            return PyLatticeDirection(self.getLatticeDirectionInC(v.lv));
        });
    }
}
#endif //OILAB_BICRYSTAL_BINDINGS_H
