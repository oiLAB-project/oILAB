//
// Created by Nikhil Chandra Admal on 7/3/25.
//

#ifndef OILAB_LATTICEDIRECTIONBINDINGS_H
#define OILAB_LATTICEDIRECTIONBINDINGS_H

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <LatticeModule.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <LatticeVectorBindings.h>




namespace pyoilab{
    template<int dim>
    class PyLatticeDirection
    {
        using PyLatticeVector = pyoilab::PyLatticeVector<dim>;

        using Lattice = gbLAB::Lattice<dim>;
        using LatticeDirection = gbLAB::LatticeDirection<dim>;
        using LatticeVector = gbLAB::LatticeVector<dim>;

        using IntScalarType = long long int;
        using MatrixDimD = Eigen::Matrix<double, dim, dim>;
        using VectorDimD = Eigen::Matrix<double, dim, 1>;
        using VectorDimI = Eigen::Matrix<IntScalarType, dim, 1>;
        using MatrixDimI = Eigen::Matrix<IntScalarType, dim, dim>;
    public:
        LatticeDirection ld;
        PyLatticeDirection(const PyLatticeVector& pv) : ld(pv.lv){}
        PyLatticeDirection(const LatticeDirection& ld) : ld(ld){}
        PyLatticeDirection(const PyLatticeDirection& pld) : ld(pld.ld){}
        PyLatticeDirection(const VectorDimI& v, const Lattice& lat) : ld(v,lat){}
        const LatticeVector& latticeVector() const
        {
            return ld.latticeVector();
        }

        VectorDimD cartesian() const {
            return ld.cartesian();
        }

        VectorDimI integerCoordinates() const {
            return ld.latticeVector();
        }

        IntScalarType dot(const PyReciprocalLatticeVector<dim>& other) const {
            return ld.dot(other.rlv);
        }
    };

    template<int dim>
    void bind_LatticeDirection(py::module_ &m) {
        using MatrixDimD = Eigen::Matrix<double, dim, dim>;
        using VectorDimD = Eigen::Matrix<double, dim, 1>;
        using VectorDimI = Eigen::Matrix<long long int, dim, 1>;

        using Lattice = gbLAB::Lattice<dim>;
        using LatticeDirection = gbLAB::LatticeDirection<dim>;
        using LatticeVector = gbLAB::LatticeVector<dim>;

        using PyLatticeDirection = PyLatticeDirection<dim>;
        using PyLatticeVector = PyLatticeVector<dim>;

        py::class_<PyLatticeDirection>(m, ("LatticeDirection" + std::to_string(dim) + "D").c_str())
                .def(py::init<const PyLatticeVector&>())
                .def(py::init<const PyLatticeDirection&>())
                .def(py::init<const VectorDimI&, const Lattice&>())
                .def("latticeVector",[](const PyLatticeDirection& ld) {
                    return PyLatticeVector(ld.latticeVector());
                },"Get the lattice vector")
                .def("cartesian",&PyLatticeDirection::cartesian, "Cartesian coordinates of the lattice direction.")
                .def("integerCoordinates",&PyLatticeDirection::integerCoordinates,"Integer coordinates of the lattice direction.")
                .def("dot",&PyLatticeDirection::dot,"dot product with a reciprocal lattice vector");
    }
}
#endif //OILAB_LATTICEDIRECTIONBINDINGS_H
