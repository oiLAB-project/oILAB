//
// Created by Nikhil Chandra Admal on 7/4/23.
//

#ifndef OILAB_PERIODICFUNCTION_H
#define OILAB_PERIODICFUNCTION_H

#include <LatticeModule.h>
#include <Eigen/Dense>
#include <MultiLattice.h>
#include <unsupported/Eigen/CXX11/Tensor>
#include <Polynomial.h>
#include <FFT.h>

using dcomplex = std::complex<double>;

namespace gbLAB {

    template<typename Scalar, int dim>
    class PeriodicFunction {
    public:
        const Lattice<dim>& lattice;
        Eigen::Tensor<Scalar,dim> values;

        // Construct a zero periodic function
        explicit PeriodicFunction(const Eigen::array<Eigen::Index,dim>& n, const Lattice<dim>& _lattice) :
                values(n), lattice(_lattice)
        {
            values.setZero();
        }

        template<typename T, int dm=dim, typename = std::enable_if_t<dm==2>>
        PeriodicFunction(const Eigen::array<Eigen::Index,dim>& n,
                         const Lattice<dim>& _lattice,
                         const Function<T,Scalar,dim>& fun) :
            values(n), lattice(_lattice)
        {
            Eigen::Vector<double,dim> center=lattice.latticeBasis.col(0)/2 + lattice.latticeBasis.col(1)/2;
            for (int i = 0; i < n[0]; i++)
            {
                for (int j = 0; j < n[1]; j++)
                {
                    Eigen::Vector<double,dim> x=i*lattice.latticeBasis.col(0)/n[0] + j*lattice.latticeBasis.col(1)/n[1];
                    values(i,j)= fun(x-center);
                }
            }
        }

        template<typename S, typename T, int dm=dim, typename = std::enable_if_t<dm==3>>
        PeriodicFunction(const Eigen::array<Eigen::Index,dim>& n,
                         const Lattice<dim>& _lattice,
                         const Function<T,Scalar,dim>& fun) :
                values(n), lattice(_lattice)
        {
            Eigen::Vector<double,dim> center=lattice.latticeBasis.col(0)/2 + lattice.latticeBasis.col(1)/2;
            for (int i = 0; i < n[0]; i++) {
                for (int j = 0; j < n[1]; j++) {
                    for (int k = 0; k < n[2]; k++) {
                        Eigen::Vector<double,dim> x= i*lattice.latticeBasis.col(0)/n[0] +
                                                     j*lattice.latticeBasis.col(1)/n[1] +
                                                     k*lattice.latticeBasis.col(2)/n[2];
                        values(i, j, k) = fun(x-center);
                    }
                }
            }
        }

        static PeriodicFunction<Scalar,dim> kernelConvolution(const Eigen::array<Eigen::Index,dim>& n,
                                                              const MultiLattice<dim>& m,
                                                              const Function<PiecewisePolynomial<Scalar,dim>,Scalar,dim>& kernel)
        {
            // create a zero periodic function
            PeriodicFunction<Scalar,dim> output(n, (const Lattice<dim>&)(m));

            // Create a periodic kernel function and compute its FFT
            const int factor=3;
            Lattice<dim> scaledLattice(factor*m.latticeBasis);
            Eigen::array<Eigen::Index,dim> nRefined;
            for(int i=0; i<dim; i++)
            {
                nRefined[i]=factor*n[i];

            }

            // pkf - periodic kernel function
            // pkf is centered w.r.t the scaled lattice
            PeriodicFunction<Scalar,dim> pkf(nRefined, scaledLattice, kernel);

            // compute the FFT of pkf
            Lattice<dim> scaledReciprocalLattice(scaledLattice.reciprocalBasis);
            PeriodicFunction<dcomplex,dim> pkfhat(nRefined,scaledReciprocalLattice);
            FFT::fft(pkf.values.template cast<dcomplex>(),pkfhat.values);

            // pkfhatCorrected - correction to pkf due to translation
            Eigen::Vector<double,dim> center=scaledLattice.latticeBasis.colwise().sum()/2;
            Exponential<dim> exp_kx(center);
            PeriodicFunction<dcomplex,dim> exponential(nRefined,scaledReciprocalLattice,exp_kx);
            PeriodicFunction<dcomplex,dim> pkfhatCorrected(nRefined,scaledReciprocalLattice);
            pkfhatCorrected.values= (exponential.values*pkfhat.values).eval();

            // Condense pkfhatCorrected
            Lattice<dim> reciprocalLattice(m.reciprocalBasis);
            PeriodicFunction<dcomplex,dim> pkfhatReduced(n,reciprocalLattice);
            Eigen::array<Eigen::DenseIndex, dim> strides;
            strides.fill(factor);
            pkfhatReduced.values= pkfhatCorrected.values.stride(strides);


            // Fourier transform of the multilattice
            PeriodicFunction<dcomplex,dim> mhat(n,reciprocalLattice);
            for(const auto& atom :  m.basisAtoms.colwise()) {
                Exponential<dim> exp_m(-atom);
                PeriodicFunction<dcomplex,dim> temp(n,reciprocalLattice,exp_m);
                mhat.values+= temp.values;
            }


            FFT::ifft(mhat.values*pkfhatReduced.values,output.values);
            return output;
        }
    };


}

#endif //OILAB_PERIODICFUNCTION_H
