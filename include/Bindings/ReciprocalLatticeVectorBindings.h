//
// Created by Nikhil Chandra Admal on 7/3/25.
//

#ifndef OILAB_RECIPROCALLATTICEVECTORBINDINGS_H
#define OILAB_RECIPROCALLATTICEVECTORBINDINGS_H

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

#include <LatticeModule.h>
#include <PyLatticeModule.h>

namespace py = pybind11;

namespace gbLAB {
// Had to wrap the ReciprocalLatticeVector<dim> class due to Pybind11 issues with classes with
// Eigen base classes
    template<int dim>
    class PyReciprocalLatticeVector {
        using Lattice = gbLAB::Lattice<dim>;
        using ReciprocalLatticeVector = gbLAB::ReciprocalLatticeVector<dim>;
        using IntScalarType = long long int;
        using MatrixDimD = Eigen::Matrix<double, dim, dim>;
        using VectorDimD = Eigen::Matrix<double, dim, 1>;
        using VectorDimI = Eigen::Matrix<IntScalarType, dim, 1>;
        using MatrixDimI = Eigen::Matrix<IntScalarType, dim, dim>;
    public:
        ReciprocalLatticeVector rlv;

        PyReciprocalLatticeVector(const Lattice &lattice) : rlv(lattice) {}

        PyReciprocalLatticeVector(const VectorDimD &cartesianCoordinates, const Lattice &lattice) : rlv(
                cartesianCoordinates,
                lattice) {}

        PyReciprocalLatticeVector(const VectorDimI &integerCoordinates, const Lattice &lattice) : rlv(
                integerCoordinates,
                lattice) {}

        PyReciprocalLatticeVector(const PyReciprocalLatticeVector &) = default;

        PyReciprocalLatticeVector(const ReciprocalLatticeVector &rlv) : rlv(rlv) {}

        PyReciprocalLatticeVector operator+(const PyReciprocalLatticeVector &other) const {
            return PyReciprocalLatticeVector(rlv.operator+(other.rlv));
        }

        PyReciprocalLatticeVector operator-(const PyReciprocalLatticeVector &other) const {
            return PyReciprocalLatticeVector(rlv.operator-(other.rlv));
        }

        PyReciprocalLatticeVector operator+=(const PyReciprocalLatticeVector &other) {
            rlv.operator+=(other.rlv);
            return *this;
        }

        PyReciprocalLatticeVector operator-=(const PyReciprocalLatticeVector &other) {
            rlv.operator-=(other.rlv);
            return *this;
        }

        PyReciprocalLatticeVector operator*(const IntScalarType &scalar) const {
            return PyReciprocalLatticeVector(rlv.operator*(scalar));
        }

        VectorDimD cartesian() const {
            return rlv.cartesian();
        }

        VectorDimI integerCoordinates() const {
            return rlv;
        }

        IntScalarType dot(const PyLatticeVector<dim> &other) const {
            return rlv.dot(other.lv);
        }

        IntScalarType planeIndexOfPoint(const VectorDimD &p) {
            return rlv.planeIndexOfPoint(p);
        }

        IntScalarType planeIndexOfPoint(const LatticeVector<dim> &p) const {
            return rlv.planeIndexOfPoint(p);
        }

        IntScalarType closestPlaneIndexOfPoint(const VectorDimD &p) const {
            return rlv.closestPlaneIndexOfPoint(p);
        }

        template<int dm=dim>
        typename std::enable_if<dm==2,PyLatticeDirection<dm>>::type
        cross(const PyReciprocalLatticeVector<dm>& other){
            return PyLatticeDirection(rlv.cross(other.rlv));
        }
        template<int dm=dim>
        typename std::enable_if<dm==3,PyLatticeDirection<dm>>::type
        cross(const PyReciprocalLatticeVector<dm>& other){
            return PyLatticeDirection(rlv.cross(other.rlv));
        }
    };

    template<int dim>
    PyReciprocalLatticeVector<dim> operator*(const long long int& scalar, const PyReciprocalLatticeVector<dim>& L) {
        return L * scalar;
    }


    template<int dim>
    void bind_ReciprocalLatticeVector(py::module_ &m) {
        using Lattice = gbLAB::Lattice<dim>;
        using LatticeVector = gbLAB::LatticeVector<dim>;
        using ReciprocalLatticeVector = gbLAB::ReciprocalLatticeVector<dim>;
        using IntScalarType = long long int;
        using MatrixDimD = Eigen::Matrix<double, dim, dim>;
        using VectorDimD = Eigen::Matrix<double, dim, 1>;
        using VectorDimI = Eigen::Matrix<IntScalarType, dim, 1>;
        using MatrixDimI = Eigen::Matrix<IntScalarType, dim, dim>;
        using PyReciprocalLatticeVector = PyReciprocalLatticeVector<dim>;

        py::class_<PyReciprocalLatticeVector>(m, ("PyReciprocalLatticeVector" + std::to_string(dim) + "D").c_str())
                .def(py::init<const Lattice &>())
                .def(py::init<const VectorDimD&, const Lattice&>())
                .def(py::init<const VectorDimI&, const Lattice&>())
                .def(py::init<const PyReciprocalLatticeVector&>())
                .def("cartesian", &PyReciprocalLatticeVector::cartesian)
                .def("integerCoordinates", &PyReciprocalLatticeVector::integerCoordinates)
                .def(py::self + py::self)
                .def(py::self - py::self)
                .def("__mul__", [](const PyReciprocalLatticeVector& self, const IntScalarType& scalar) {
                    return self * scalar;
                }, py::is_operator())
                .def("__rmul__", [](const PyReciprocalLatticeVector& self, const IntScalarType& scalar) {
                    return scalar * self;
                }, py::is_operator())
                .def(py::self -= py::self)
                .def(py::self += py::self)
                .def("dot",&PyReciprocalLatticeVector::dot)
                // note that cross is a template member function
                .def("cross",&PyReciprocalLatticeVector::template cross<dim>);
    }

}
#endif //OILAB_RECIPROCALLATTICEVECTORBINDINGS_H