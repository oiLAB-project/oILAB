//
// Created by Nikhil Chandra Admal on 5/31/24.
//
//#include <LatticeFunction.h>
#include <FFT.h>

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
    PeriodicFunction<typename LatticeFunction<Scalar,dim>::dcomplex,dim> LatticeFunction<Scalar,dim>::ifft()
    {
        Eigen::Matrix<double,Eigen::Dynamic,dim> unitCell(basisVectors.transpose().completeOrthogonalDecomposition().pseudoInverse());
        PeriodicFunction<dcomplex,dim> pf(values.dimensions(),unitCell);
        FFT::ifft(values.template cast<dcomplex>(),pf.values);
        return pf;
    }

}
