#include <pybind11/pybind11.h>
#include <LatticeModule.h>
#include <Eigen/Dense>
#include <pybind11/eigen.h>

namespace py = pybind11;

void bind_Lattice(py::module_& m)
{
    using Lattice2D = gbLAB::Lattice<2>;
    using ReciprocalLatticeDirection2D = gbLAB::ReciprocalLatticeDirection<2>;
    using Lattice3D = gbLAB::Lattice<3>;
    using ReciprocalLatticeDirection3D = gbLAB::ReciprocalLatticeDirection<3>;

    py::class_<Lattice2D>(m, "Lattice2D")
        .def(py::init<const Eigen::Matrix2d&, Eigen::Matrix2d&>(), py::arg("A"), py::arg("Q") = Eigen::Matrix2d::Identity())
        .def("interPlanarSpacing", &Lattice2D::interPlanarSpacing);

    py::class_<Lattice3D>(m, "Lattice3D")
        .def(py::init<const Eigen::Matrix3d&, const Eigen::Matrix3d&>(), py::arg("A"), py::arg("Q") = Eigen::Matrix3d::Identity())
        .def("interPlanarSpacing", &Lattice3D::interPlanarSpacing);
}

