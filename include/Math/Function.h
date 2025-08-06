//
// Created by Nikhil Chandra Admal on 7/9/23.
//

#ifndef OILAB_POLYNOMIAL_H
#define OILAB_POLYNOMIAL_H

#include <deque>
#include <cmath>
#include "Eigen/Dense"
#include "unsupported/Eigen/CXX11/Tensor"
#include <array>

namespace gbLAB {

    template<typename Scalar, int dim>
    class LatticeFunction;

    template<typename Scalar, int dim>
    class PeriodicFunction;

    template<typename Derived, typename Scalar>
    class Function {
        using dcomplex= std::complex<double>;
    private:
        const Derived& derivedFunction;
    public:
        double domainSize;
        explicit Function(double _domainSize = std::numeric_limits<double>::infinity());
        Scalar operator()(const Eigen::Vector<double,Eigen::Dynamic>& vec) const;

        /*
        template<int dim>
        LatticeFunction<dcomplex,dim> fft(const std::array<Eigen::Index,dim>& n, const Eigen::Matrix<double,Eigen::Dynamic,dim>& basisVectors) const;
         */
    };

    /* ******************************************** */
    class Exponential : public Function<Exponential,std::complex<double>>
    {
    private:
        const Eigen::Vector<double,Eigen::Dynamic>& x;
    public:
        explicit Exponential(const Eigen::Vector<double,Eigen::Dynamic>& _x);
        std::complex<double> operator()(const Eigen::Vector<double,Eigen::Dynamic>& vec) const;
    };

    /* ******************************************** */
    // it is a 2d scalar function that takes 3d input
    template<typename T, typename Scalar>
    class Shift : public Function<Shift<T,Scalar>, Scalar>
    {
    public:
        const Function<T,Scalar>& fun;
        Eigen::Vector<double,Eigen::Dynamic> t;
        Shift(const Eigen::Vector<double,Eigen::Dynamic>& t, const Function<T,Scalar>& fun);
        Scalar operator()(const Eigen::Vector<double,Eigen::Dynamic>& y) const;
    };
}

#include <FunctionImplementation.h>
#endif //OILAB_POLYNOMIAL_H
