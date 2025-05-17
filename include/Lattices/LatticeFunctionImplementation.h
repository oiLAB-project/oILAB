//
// Created by Nikhil Chandra Admal on 5/31/24.
//
#ifndef OILAB_LATTICEFUNCTIONIMPLEMENTATION_H
#define OILAB_LATTICEFUNCTIONIMPLEMENTATION_H
#include <FFT.h>
#include <iostream>

namespace gbLAB
{
    template<typename Scalar, int dim>
    LatticeFunction<Scalar,dim>::LatticeFunction(const Eigen::array<Eigen::Index,dim>& n,
                                                 const Eigen::Matrix<double,Eigen::Dynamic,dim>& _basisVectors) :
            values(n), basisVectors(_basisVectors)
    {
        values.setZero();
    }

    template<typename Scalar, int dim>
    template<typename T, typename, typename, int dm, typename>
    LatticeFunction<Scalar,dim>::LatticeFunction(const Eigen::array<Eigen::Index, dim> &n,
                    const Eigen::Matrix<double, Eigen::Dynamic, dim> &_basisVectors,
                    const Function<T,Scalar>& fun) :
            values(n), basisVectors(_basisVectors) {
        for (int i = 0; i < n[0]; i++) {
            int in= i > n[0]/2 ? i-n[0] : i;
            values(i) = fun(in * basisVectors.col(0));
        }
    }

    template<typename Scalar, int dim>
    template<typename T, typename, int dm, typename>
    LatticeFunction<Scalar,dim>::LatticeFunction(const Eigen::array<Eigen::Index, dim> &n,
                                                 const Eigen::Matrix<double,Eigen::Dynamic,dim>& _basisVectors,
                                                 const Function<T,Scalar>& fun) :
            values(n), basisVectors(_basisVectors) {
        for (int i = 0; i < n[0]; i++) {
            for (int j = 0; j < n[1]; j++) {
                //values(i, j) = fun(i * basisVectors.col(0) + j * basisVectors.col(1));
                int in= i > n[0]/2 ? i-n[0] : i;
                int jn= j > n[1]/2 ? j-n[1] : j;
                values(i,j) = fun(in * basisVectors.col(0) + jn * basisVectors.col(1));
            }
        }
    }

    template<typename Scalar, int dim>
    template<typename T, int dm, typename>
    LatticeFunction<Scalar,dim>::LatticeFunction(const Eigen::array<Eigen::Index, dim> &n,
                    const Eigen::Matrix<double, Eigen::Dynamic, dim> &_basisVectors,
                    const Function<T, Scalar> &fun) :
            values(n), basisVectors(_basisVectors) {
        for (int i = 0; i < n[0]; i++) {
            for (int j = 0; j < n[1]; j++) {
                for (int k = 0; k < n[2]; k++) {
                    //values(i, j, k) = fun(i * basisVectors.col(0) + j * basisVectors.col(1) + k * basisVectors.col(2));
                    int in= i>n[0]/2 ? i-n[0] : i;
                    int jn= j>n[1]/2 ? j-n[1] : j;
                    int kn= k>n[2]/2 ? k-n[2] : k;
                    values(i, j, k) = fun(in * basisVectors.col(0) +
                                          jn * basisVectors.col(1) +
                                          kn * basisVectors.col(2));
                }
            }
        }
    }

    template<typename Scalar, int dim>
    std::complex<double> LatticeFunction<Scalar,dim>::dot(const LatticeFunction<std::complex<double>,dim>& other) const
    {
        Eigen::Tensor<std::complex<double>,0> sum((this->values * other.values.conjugate()).sum());
        Eigen::Matrix<double,dim,dim> gramMatrix;
        for(int i=0; i<dim; ++i)
            for(int j=0; j<dim; ++j)
                gramMatrix(i,j)= basisVectors.col(i).dot(basisVectors.col(j));
        return sum(0) * sqrt(gramMatrix.determinant());

    }

    template<typename Scalar, int dim>
    PeriodicFunction<typename LatticeFunction<Scalar,dim>::dcomplex,dim> LatticeFunction<Scalar,dim>::ifft() const
    {
        Eigen::Matrix<double,Eigen::Dynamic,dim> unitCell(basisVectors.transpose().completeOrthogonalDecomposition().pseudoInverse());
        PeriodicFunction<dcomplex,dim> pf(values.dimensions(),unitCell);
        FFT::ifft(values.template cast<dcomplex>(),pf.values);
        //return pf;
        // Calculate the area spanned by the unit cell vectors
        Eigen::Matrix<double,dim,dim> unitCellGramMatrix;
        for(int i=0; i<dim; ++i)
            for(int j=0; j<dim; ++j)
                unitCellGramMatrix(i,j)= unitCell.col(i).dot(unitCell.col(j));

        Eigen::array<Eigen::Index,dim> n= this->values.dimensions();
        int prod= std::accumulate(std::begin(n),
                                  std::begin(n) + dim,
                                  1,
                                  std::multiplies<>{});
        pf.values=  pf.values * (dcomplex)(prod/sqrt(unitCellGramMatrix.determinant()));
        return pf;
    }

    template<typename Scalar, int dim>
    LatticeFunction<Scalar, dim> operator*(const LatticeFunction<Scalar,dim>& lf1, const LatticeFunction<Scalar,dim>& lf2)
    {
        // assert that the two lattice functions have the same domain
        const auto n= lf1.values.dimensions();
        LatticeFunction<Scalar,dim> output(n,lf1.basisVectors);
        output.values= lf1.values * lf2.values;
        return output;
    }
}
#endif // OILAB_LATTICEFUNCTIONIMPLEMENTATION_H