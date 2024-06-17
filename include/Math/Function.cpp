//
// Created by Nikhil Chandra Admal on 5/30/24.
//

//#include <Function.h>
//#include <PeriodicFunction.h>
#include <iostream>

namespace gbLAB
{
    template<typename Derived, typename Scalar>
    Function<Derived,Scalar>::Function(double _domainSize) :
        derivedFunction( static_cast<const Derived&>(*this)),
        domainSize(_domainSize)
    {}

    template<typename Derived, typename Scalar>
    Scalar Function<Derived,Scalar>::operator()(const Eigen::Vector<double,Eigen::Dynamic>& vec) const
    {
        return derivedFunction(vec);
    }

    template<typename Derived, typename Scalar>
    template<int dim>
    LatticeFunction<typename Function<Derived,Scalar>::dcomplex,dim>
            Function<Derived,Scalar>::fft(const std::array<Eigen::Index,dim>& n, const Eigen::Matrix<double,Eigen::Dynamic,dim>& basisVectors) const
    {
        // Initialize a zero periodic function and its fft
        LatticeFunction<dcomplex,dim> pfhat(n,basisVectors);
        PeriodicFunction<std::complex<double>,dim> pf_complex(pfhat.ifft());
        PeriodicFunction<double,dim> pf(pf_complex.values.dimensions(),pf_complex.unitCell);
        pf.values= pf_complex.values.real();

        Eigen::Vector<double, Eigen::Dynamic> factor, inversePlaneDistances;
        inversePlaneDistances = basisVectors.colwise().norm();
        factor = (2.0 * domainSize * inversePlaneDistances).array().ceil();
        //std::cout << "Scaling the lattice by factors = " << factor.transpose() << std::endl;
        Eigen::Matrix<double, Eigen::Dynamic, dim> scaledUnitCell= pf.unitCell.array().rowwise() * factor.array().transpose();

        Eigen::array<Eigen::Index, dim> nRefined = n;
        for (int i = 0; i < dim; i++) {
            nRefined[i] = factor(i) * n[i];
        }

        // pkf - periodic kernel function
        PeriodicFunction<double,dim> pfRefined(nRefined, scaledUnitCell, *this);
        LatticeFunction<dcomplex, dim> pfhatRefined(pfRefined.fft());

        /*
        // Correct pkfkhat as the kernel function was centered w.r.t the scaled lattice
        LatticeFunction<dcomplex,dim> temp(nRefined, pfhatRefined.basisVectors);
        Eigen::Vector<double, Eigen::Dynamic> center = pfRefined.unitCell.rowwise().sum()/2;
        Exponential exp_function(center);
        LatticeFunction<dcomplex, dim> exp_factor(nRefined, pfhatRefined.basisVectors, exp_function);
        temp.values = exp_factor.values * pfhatRefined.values;
        pfhatRefined.values = temp.values;
        */

        Eigen::array<Eigen::DenseIndex, dim> strides;
        strides.fill(1);
        for (int i = 0; i < dim; i++)
            strides[i] = factor(i);

        pfhat.values = pfhatRefined.values.stride(strides);

        return pfhat;

    }

    /* ******************************************** */

    Exponential::Exponential(const Eigen::Vector<double,Eigen::Dynamic>& _x) : x(_x){}
    std::complex<double> Exponential::operator()(const Eigen::Vector<double,Eigen::Dynamic>& vec) const
    {
        return exp(2*M_PI*std::complex<double>(0,1)*vec.dot(x));
    }

    /* ******************************************** */
    template<typename T, typename Scalar>
    Shift<T,Scalar>::Shift(const Eigen::Vector<double,Eigen::Dynamic>& t,
                           const Function<T,Scalar>& fun):
            Function<Shift<T,Scalar>,Scalar>(fun.domainSize),
            t(t),fun(fun)
    {}
    template<typename T, typename Scalar>
    Scalar Shift<T,Scalar>::operator()(const Eigen::Vector<double,Eigen::Dynamic>& y) const
    {
        //return fun(y+t);
        return fun(t-y);
    }
}
