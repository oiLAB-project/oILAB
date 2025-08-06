//
// Created by Nikhil Chandra Admal on 7/4/25.
//

#ifndef OILAB_RECIPROCALLATTICEDIRECTIONBINDINGS_H
#define OILAB_RECIPROCALLATTICEDIRECTIONBINDINGS_H

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <LatticeModule.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

namespace pyoilab{
    template<int dim>
    class PyReciprocalLatticeDirection
    {
        using Lattice = gbLAB::Lattice<dim>;
        using ReciprocalLatticeDirection = gbLAB::ReciprocalLatticeDirection<dim>;
        using ReciprocalLatticeVector = gbLAB::ReciprocalLatticeVector<dim>;
        using PyReciprocalLatticeVector = pyoilab::PyReciprocalLatticeVector<dim>;

        using IntScalarType = long long int;
        using MatrixDimD = Eigen::Matrix<double, dim, dim>;
        using VectorDimD = Eigen::Matrix<double, dim, 1>;
        using VectorDimI = Eigen::Matrix<IntScalarType, dim, 1>;
        using MatrixDimI = Eigen::Matrix<IntScalarType, dim, dim>;
    public:
        ReciprocalLatticeDirection rld;

        PyReciprocalLatticeDirection(const PyReciprocalLatticeVector& prlv) : rld(prlv.rlv){}
        PyReciprocalLatticeDirection(const ReciprocalLatticeVector& rlv) : rld(rlv){}
        PyReciprocalLatticeDirection(const PyReciprocalLatticeDirection& prld) = default;
        PyReciprocalLatticeDirection(const ReciprocalLatticeDirection& rld) : rld(rld){}

        const ReciprocalLatticeVector& reciprocalLatticeVector() const
        {
            return rld.reciprocalLatticeVector();
        }

        VectorDimD cartesian() const {
            return rld.cartesian();
        }

        VectorDimI integerCoordinates() const {
            return rld.reciprocalLatticeVector();
        }

        IntScalarType dot(const PyLatticeVector<dim>& other) const {
            return rld.dot(other.lv);
        }

        double planeSpacing() const {
            return rld.planeSpacing();
        }

        int stacking() const {
            return rld.stacking();
        }
    };

    template<int dim>
    void bind_ReciprocalLatticeDirection(py::module_ &m) {
        using MatrixDimD = Eigen::Matrix<double, dim, dim>;
        using VectorDimD = Eigen::Matrix<double, dim, 1>;
        using VectorDimI = Eigen::Matrix<long long int, dim, 1>;

        using Lattice = gbLAB::Lattice<dim>;
        using ReciprocalLatticeDirection = gbLAB::ReciprocalLatticeDirection<dim>;
        using LatticeVector = gbLAB::LatticeVector<dim>;

        using PyReciprocalLatticeDirection = PyReciprocalLatticeDirection<dim>;
        using PyReciprocalLatticeVector = PyReciprocalLatticeVector<dim>;

        py::class_<PyReciprocalLatticeDirection>(m, ("ReciprocalLatticeDirection" + std::to_string(dim) + "D").c_str())
                .def(py::init<const PyReciprocalLatticeVector&>())
                .def(py::init<const PyReciprocalLatticeDirection&>())
                .def("reciprocalLatticeVector",[](const PyReciprocalLatticeDirection& rld) {
                    return PyReciprocalLatticeVector(rld.reciprocalLatticeVector());
                },"Get the reciprocal lattice vector")
                .def("cartesian",&PyReciprocalLatticeDirection::cartesian, "Cartesian coordinates of the reciprocal lattice direction.")
                .def("integerCoordinates",&PyReciprocalLatticeDirection::integerCoordinates,"Integer coordinates of the reciprocal lattice direction.")
                .def("dot",&PyReciprocalLatticeDirection::dot,"dot product with a lattice vector")
                .def("planeSpacing",&PyReciprocalLatticeDirection::planeSpacing,"distance between two consecutive planes")
                .def("stacking",&PyReciprocalLatticeDirection::stacking,"stacking");

    }
}
#endif //OILAB_RECIPROCALLATTICEDIRECTIONBINDINGS_H
