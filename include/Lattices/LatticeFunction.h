//
// Created by Nikhil Chandra Admal on 5/26/24.
//

#ifndef OILAB_LATTICEFUNCTION_H
#define OILAB_LATTICEFUNCTION_H

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

// Lattice function defined on the lattice points
namespace gbLAB {

    template<typename Derived, typename Scalar>
    class Function;

    template<typename Scalar, int dim>
    class PeriodicFunction;

    template<typename Scalar, int dim>
    class LatticeFunction {
        using dcomplex= std::complex<double>;
    public:
        const Eigen::Matrix<double, Eigen::Dynamic, dim> basisVectors;
        Eigen::Tensor<Scalar, dim> values;
        explicit LatticeFunction(const Eigen::array<Eigen::Index, dim>& n,
                                 const Eigen::Matrix<double, Eigen::Dynamic, dim>& _basisVectors);

        template<typename T, typename=T, typename=T, int dm = dim, typename = std::enable_if_t<dm == 1>>
        LatticeFunction(const Eigen::array<Eigen::Index, dim> &n,
                        const Eigen::Matrix<double, Eigen::Dynamic, dim> &_basisVectors,
                        const Function<T,Scalar>& fun);

        template<typename T, typename=T, int dm = dim, typename = std::enable_if_t<dm == 2>>
        LatticeFunction(const Eigen::array<Eigen::Index, dim> &n,
                        const Eigen::Matrix<double, Eigen::Dynamic, dim> &_basisVectors,
                        const Function<T, Scalar> &fun);

        template<typename T, int dm = dim, typename = std::enable_if_t<dm == 3>>
        LatticeFunction(const Eigen::array<Eigen::Index, dim> &n,
                        const Eigen::Matrix<double, Eigen::Dynamic, dim> &_basisVectors,
                        const Function<T, Scalar> &fun);

        std::complex<double> dot(const LatticeFunction<std::complex<double>,dim>& other) const;

        PeriodicFunction<dcomplex,dim> ifft() const;
    };

    template<typename Scalar, int dim>
    LatticeFunction<Scalar, dim> operator*(const LatticeFunction<Scalar,dim>& lf1, const LatticeFunction<Scalar,dim>& lf2);
}

#include <LatticeFunctionImplementation.h>

#endif //OILAB_LATTICEFUNCTION_H
