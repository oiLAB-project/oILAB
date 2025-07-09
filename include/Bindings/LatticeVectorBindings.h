//
// Created by Nikhil Chandra Admal on 7/3/25.
//

#ifndef OILAB_LATTICEVECTORBINDINGS_H
#define OILAB_LATTICEVECTORBINDINGS_H

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <LatticeModule.h>
#include <PyLatticeModule.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

namespace pyoilab{
// Had to wrap the LatticeVector<dim> class due to Pybind11 issues with classes with
// Eigen base classes
    template<int dim>
    class PyLatticeVector {
        using Lattice = gbLAB::Lattice<dim>;
        using LatticeVector = gbLAB::LatticeVector<dim>;
        using PyReciprocalLatticeDirection = PyReciprocalLatticeDirection<dim>;


        using IntScalarType = long long int;
        using MatrixDimD = Eigen::Matrix<double, dim, dim>;
        using VectorDimD = Eigen::Matrix<double, dim, 1>;
        using VectorDimI = Eigen::Matrix<IntScalarType, dim, 1>;
        using MatrixDimI = Eigen::Matrix<IntScalarType, dim, dim>;
    public:
        LatticeVector lv;

        PyLatticeVector(const Lattice& lattice) : lv(lattice) {}

        PyLatticeVector(const VectorDimD& cartesianCoordinates, const Lattice& lattice) : lv(cartesianCoordinates,
                                                                                             lattice) {}

        PyLatticeVector(const VectorDimI& integerCoordinates, const Lattice& lattice) : lv(integerCoordinates,
                                                                                           lattice) {}

        PyLatticeVector(const PyLatticeVector&) = default;

        PyLatticeVector(const LatticeVector& lv) : lv(lv) {}

        PyLatticeVector operator+(const PyLatticeVector& other) const {
            return PyLatticeVector(lv.operator+(other.lv));
        }

        PyLatticeVector operator-(const PyLatticeVector& other) const {
            return PyLatticeVector(lv.operator-(other.lv));
        }

        PyLatticeVector operator+=(const PyLatticeVector& other) {
            lv.operator+=(other.lv);
            return *this;
        }

        PyLatticeVector operator-=(const PyLatticeVector& other) {
            lv.operator-=(other.lv);
            return *this;
        }

        PyLatticeVector operator*(const IntScalarType& scalar) const {
            return PyLatticeVector(lv.operator*(scalar));
        }

        VectorDimD cartesian() const {
            return lv.cartesian();
        }

        VectorDimI integerCoordinates() const {
            return lv;
        }

        void integerCoordinates(const VectorDimI& input) {
            lv << input;
        }

        IntScalarType dot(const PyReciprocalLatticeVector<dim>& other){
            return lv.dot(other.rlv);
        }

        template<int dm=dim>
        typename std::enable_if<dm==2,PyReciprocalLatticeDirection>::type
        cross(const PyLatticeVector<dim>& other) const
        {
            return PyReciprocalLatticeDirection(lv.cross(other.lv));
        }

        template<int dm=dim>
        typename std::enable_if<dm==3,PyReciprocalLatticeDirection>::type
        cross(const PyLatticeVector<dm>& other) const
        {
            return PyReciprocalLatticeDirection(lv.cross(other.lv));
        }

        template<int dm=dim>
        typename std::enable_if<dm==2,PyReciprocalLatticeDirection>::type
        cross() const
        {
            return PyReciprocalLatticeDirection(lv.cross());
        }
        template<int dm=dim>
        typename std::enable_if<dm==3,PyReciprocalLatticeDirection>::type
        cross() const
        {
            return PyReciprocalLatticeDirection(lv.cross());
        }

    };

    template<int dim>
    PyLatticeVector<dim> operator*(const long long int& scalar, const PyLatticeVector<dim>& L) {
        return L * scalar;
    }

    template<int dim>
    void bind_LatticeVector(py::module_& m) {
        using Lattice = gbLAB::Lattice<dim>;
        using LatticeVector = gbLAB::LatticeVector<dim>;
        using IntScalarType = long long int;
        using MatrixDimD = Eigen::Matrix<double, dim, dim>;
        using VectorDimD = Eigen::Matrix<double, dim, 1>;
        using VectorDimI = Eigen::Matrix<IntScalarType, dim, 1>;
        using MatrixDimI = Eigen::Matrix<IntScalarType, dim, dim>;
        using PyLatticeVector = PyLatticeVector<dim>;
        using PyReciprocalLatticeDirection = PyReciprocalLatticeDirection<dim>;

        py::class_<PyLatticeVector>(m, ("LatticeVector" + std::to_string(dim) + "D").c_str())
                .def(py::init<const Lattice&>())
                .def(py::init<const VectorDimD&, const Lattice&>())
                .def(py::init<const VectorDimI&, const Lattice&>())
                .def(py::init<const PyLatticeVector&>())

                .def("cartesian", &PyLatticeVector::cartesian)
                .def("integerCoordinates", static_cast<VectorDimI (PyLatticeVector::*)() const>(&PyLatticeVector::integerCoordinates),
                     "output the integer coordinates of the lattice vecctor")
                .def("integerCoordinates", static_cast<void (PyLatticeVector::*)(const VectorDimI&)>(&PyLatticeVector::integerCoordinates),
                    "input the integer coordinates of the lattice vector")
                .def("dot", &PyLatticeVector::dot,
                     "dot product with a reciprocal lattice vector")
                .def(py::self + py::self)
                .def(py::self - py::self)
                .def("__mul__", [](const PyLatticeVector& self, const IntScalarType& scalar) {
                    return self * scalar;
                }, py::is_operator())
                .def("__rmul__", [](const PyLatticeVector& self, const IntScalarType& scalar) {
                    return scalar * self;
                }, py::is_operator())
                .def(py::self -= py::self)
                .def(py::self += py::self)
                .def("cross", static_cast<PyReciprocalLatticeDirection (PyLatticeVector::*)(const PyLatticeVector&) const> (&PyLatticeVector::cross))
                .def("cross", static_cast<PyReciprocalLatticeDirection (PyLatticeVector::*)() const> (&PyLatticeVector::cross));
    }
}
#endif //OILAB_LATTICEVECTORBINDINGS_H
