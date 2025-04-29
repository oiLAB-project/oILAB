//
// Created by Nikhil Chandra Admal on 5/31/24.
//
#ifndef OILAB_PERIODICFUNCTIONIMPLEMENTATION_H
#define OILAB_PERIODICFUNCTIONIMPLEMENTATION_H
//#include <PeriodicFunction.h>
#include <FFT.h>

namespace gbLAB
{
    template<typename Scalar, int dim>
    PeriodicFunction<Scalar,dim>::PeriodicFunction(const Eigen::array<Eigen::Index,dim>& n,
                                                   const Eigen::Matrix<double,Eigen::Dynamic,dim>&  _unitCell) :
            values(n), unitCell(_unitCell)
    {
        values.setZero();
    }

    template<typename Scalar, int dim>
    template<typename T, typename, typename, int dm, typename>
    PeriodicFunction<Scalar,dim>::PeriodicFunction(const Eigen::array<Eigen::Index,dim>& n,
                     const Eigen::Matrix<double,Eigen::Dynamic,dim>&  _unitCell,
                     const Function<T,Scalar>& fun) :
            values(n), unitCell(_unitCell)
    {
        //Eigen::Vector<double,Eigen::Dynamic> center= unitCell.col(0)/2;
        Eigen::Vector<double,Eigen::Dynamic> center= unitCell.col(0)/2;
        for (int i = 0; i < n[0]; i++)
        {
            Eigen::Vector<double,Eigen::Dynamic> x=i*unitCell.col(0)/n[0];
            //values(i)= fun(x-center);
            values(i)= fun(x);
        }
    }

    template<typename Scalar, int dim>
    template<typename T, typename, int dm, typename>
    PeriodicFunction<Scalar,dim>::PeriodicFunction(const Eigen::array<Eigen::Index,dim>& n,
                     const Eigen::Matrix<double,Eigen::Dynamic,dim>&  _unitCell,
                     const Function<T,Scalar>& fun) :
            values(n), unitCell(_unitCell)
    {
        //Eigen::Vector<double,Eigen::Dynamic> center= unitCell.col(0)/2 + unitCell.col(1)/2;
        for (int i = 0; i < n[0]; i++)
        {
            for (int j = 0; j < n[1]; j++)
            {
                Eigen::Vector<double,Eigen::Dynamic> x=i*unitCell.col(0)/n[0] + j*unitCell.col(1)/n[1];
                //values(i,j)= fun(x-center);
                values(i,j)= fun(x);
            }
        }
    }

    template<typename Scalar, int dim>
    template<typename T, int dm, typename>
    PeriodicFunction<Scalar,dim>::PeriodicFunction(const Eigen::array<Eigen::Index,dim>& n,
                     const Eigen::Matrix<double,Eigen::Dynamic,dim>&  _unitCell,
                     const Function<T,Scalar>& fun) :
            values(n), unitCell(_unitCell)
    {
        //Eigen::Vector<double,Eigen::Dynamic> center=unitCell.col(0)/2 + unitCell.col(1)/2 + unitCell.col(2)/2;
        for (int i = 0; i < n[0]; i++) {
            for (int j = 0; j < n[1]; j++) {
                for (int k = 0; k < n[2]; k++) {
                    Eigen::Vector<double,Eigen::Dynamic> x= i*unitCell.col(0)/n[0] +
                                                 j*unitCell.col(1)/n[1] +
                                                 k*unitCell.col(2)/n[2];
                    //values(i, j, k) = fun(x-center);
                    values(i, j, k) = fun(x);
                }
            }
        }
    }



    template<typename Scalar, int dim>
    LatticeFunction<typename PeriodicFunction<Scalar,dim>::dcomplex,dim> PeriodicFunction<Scalar,dim>::fft() const
    {
        Eigen::Matrix<double,Eigen::Dynamic,dim> basisVectors(unitCell.transpose().completeOrthogonalDecomposition().pseudoInverse());
        LatticeFunction<dcomplex,dim> pfhat(values.dimensions(),basisVectors);
        FFT::fft(values.template cast<dcomplex>(),pfhat.values);
        return pfhat;
    }

    template<typename Scalar, int dim>
    double PeriodicFunction<Scalar,dim>::dot(const PeriodicFunction<Scalar,dim>& other) const
    {
        Eigen::Tensor<double,0> sum((this->values * other.values).sum());
        Eigen::Matrix<double,dim,dim> gramMatrix;
        for(int i=0; i<dim; ++i)
            for(int j=0; j<dim; ++j)
                gramMatrix(i,j)= unitCell.col(i).dot(unitCell.col(j));
        Eigen::array<Eigen::Index,dim> n= this->values.dimensions();
        int prod= std::accumulate(std::begin(n),
                        std::begin(n) + dim,
                        1,
                        std::multiplies<>{});
        return sum(0)*sqrt(gramMatrix.determinant())/prod;
    }

    template<typename Scalar, int dim> template <typename T>
    PeriodicFunction<Scalar,dim> PeriodicFunction<Scalar,dim>::kernelConvolution(const Function<T,Scalar>& kernel)
    {
        // pfhat - fft of the periodic function
        const auto pfhat(fft());
        Eigen::array<Eigen::Index,dim> n= values.dimensions();
        PeriodicFunction<Scalar,dim> output(n, unitCell);

        // fourier transform of the kernel function
        LatticeFunction<Scalar, dim> pkfhat(kernel.fft(n,pfhat.basisVectors));

        PeriodicFunction<dcomplex,dim> tempOutput(n,unitCell);
        FFT::ifft(pfhat.values*pkfhat.values,tempOutput.values);
        output.values = tempOutput.values.real();
        return output;
    }

}
#endif // OILAB_PERIODICFUNCTIONIMPLEMENTATION_H
