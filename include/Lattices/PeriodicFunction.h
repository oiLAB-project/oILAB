//
// Created by Nikhil Chandra Admal on 7/4/23.
//

#ifndef OILAB_PERIODICFUNCTION_H
#define OILAB_PERIODICFUNCTION_H

#include <LatticeModule.h>
#include <Eigen/Dense>
#include <MultiLattice.h>
#include <unsupported/Eigen/CXX11/Tensor>
#include <Function.h>
#include <FFT.h>

using dcomplex = std::complex<double>;

namespace gbLAB {

    template<typename Scalar, int dim>
    class LatticeFunction {
    public:
        const Lattice<dim>& lattice;
        Eigen::Tensor<Scalar,dim> values;
        // Construct a zero Lattice function
        explicit LatticeFunction(const Eigen::array<Eigen::Index,dim>& n, const Lattice<dim>& _lattice) :
                values(n), lattice(_lattice)
        {
            values.setZero();
        }

        template<typename T, typename=T,typename=T, int dm=dim, typename = std::enable_if_t<dm==1>>
        LatticeFunction(const Eigen::array<Eigen::Index,dim>& n,
                        const Lattice<dim>& _lattice,
                        const Function<T,Scalar,dim>& fun) :
                values(n), lattice(_lattice)
        {
            for (int i = 0; i < n[0]; i++)
            {
                LatticeVector<dm> x(lattice);
                x << i;
                values(i)= fun(x.cartesian());
            }
        }

        template<typename T, typename=T, int dm=dim, typename = std::enable_if_t<dm==2>>
        LatticeFunction(const Eigen::array<Eigen::Index,dim>& n,
                         const Lattice<dim>& _lattice,
                         const Function<T,Scalar,dim>& fun) :
                values(n), lattice(_lattice)
        {
            for (int i = 0; i < n[0]; i++)
            {
                for (int j = 0; j < n[1]; j++)
                {
                    LatticeVector<dm> x(lattice);
                    x << i,j;
                    values(i,j)= fun(x.cartesian());
                }
            }
        }

        template<typename T, int dm=dim, typename = std::enable_if_t<dm==3>>
        LatticeFunction(const Eigen::array<Eigen::Index,dim>& n,
                         const Lattice<dim>& _lattice,
                         const Function<T,Scalar,dim>& fun) :
                values(n), lattice(_lattice)
        {
            for (int i = 0; i < n[0]; i++) {
                for (int j = 0; j < n[1]; j++) {
                    for (int k = 0; k < n[2]; k++) {
                        LatticeVector<dm> x(lattice);
                        x << i,j,k;
                        values(i, j, k) = fun(x.cartesian());
                    }
                }
            }
        }
    };

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

        template<typename T, typename = T, typename=T, int dm=dim, typename = std::enable_if_t<dm==1>>
        //template<typename T, int dm=dim, typename = std::enable_if_t<dm==2>>
        PeriodicFunction(const Eigen::array<Eigen::Index,dim>& n,
                         const Lattice<dim>& _lattice,
                         const Function<T,Scalar,dim>& fun) :
                values(n), lattice(_lattice)
        {
            Eigen::Vector<double,dim> center= lattice.latticeBasis.col(0)/2;
            for (int i = 0; i < n[0]; i++)
            {
                Eigen::Vector<double,dim> x=i*lattice.latticeBasis.col(0)/n[0];
                values(i)= fun(x-center);
            }
        }

        template<typename T, typename = T, int dm=dim, typename = std::enable_if_t<dm==2>>
        //template<typename T, int dm=dim, typename = std::enable_if_t<dm==2>>
        PeriodicFunction(const Eigen::array<Eigen::Index,dim>& n,
                         const Lattice<dim>& _lattice,
                         const Function<T,Scalar,dim>& fun) :
            values(n), lattice(_lattice)
        {
            Eigen::Vector<double,dim> center= lattice.latticeBasis.col(0)/2 + lattice.latticeBasis.col(1)/2;
            for (int i = 0; i < n[0]; i++)
            {
                for (int j = 0; j < n[1]; j++)
                {
                    Eigen::Vector<double,dim> x=i*lattice.latticeBasis.col(0)/n[0] + j*lattice.latticeBasis.col(1)/n[1];
                    values(i,j)= fun(x-center);
                }
            }
        }

        template<typename T, int dm=dim, typename = std::enable_if_t<dm==3>>
        PeriodicFunction(const Eigen::array<Eigen::Index,dim>& n,
                         const Lattice<dim>& _lattice,
                         const Function<T,Scalar,dim>& fun) :
                values(n), lattice(_lattice)
        {
            //Eigen::Vector<double,dim> center=lattice.latticeBasis.col(0)/2 + lattice.latticeBasis.col(1)/2;
            Eigen::Vector<double,dim> center=lattice.latticeBasis.col(0)/2 +
                                             lattice.latticeBasis.col(1)/2 +
                    lattice.latticeBasis.col(2)/2;
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

        template <typename T>
        static PeriodicFunction<Scalar,dim> kernelConvolution(const Eigen::array<Eigen::Index,dim>& n,
                                                              const MultiLattice<dim>& m,
                                                              const Function<T,Scalar,dim>& kernel)
        {
            // create a zero periodic function
            PeriodicFunction<Scalar,dim> output(n, (const Lattice<dim>&)(m));
            PeriodicFunction<dcomplex,dim> tempOutput(n, (const Lattice<dim>&)(m));

            // Create a periodic kernel function and compute its FFT
            Eigen::Vector<double, dim> factor, inversePlaneDistances;
            inversePlaneDistances= m.reciprocalBasis.colwise().norm();
            factor= (2.0*kernel.domainSize * inversePlaneDistances).array().ceil();

            std::cout << "Scaling the lattice by factors = " << factor.transpose() << std::endl;
            Eigen::Matrix<double,dim,dim> basisVectors= m.latticeBasis.array().colwise() * factor.array();
            Lattice<dim> scaledLattice(basisVectors);

            Eigen::array<Eigen::Index,dim> nRefined= n;
            for(int i=0; i<dim; i++)
            {
                nRefined[i]=factor(i)*n[i];
            }

            // pkf - periodic kernel function
            // pkf is centered w.r.t the scaled lattice
            PeriodicFunction<Scalar,dim> pkf(nRefined, scaledLattice, kernel);

            // compute the FFT of pkf
            Lattice<dim> scaledReciprocalLattice(scaledLattice.reciprocalBasis);
            LatticeFunction<dcomplex,dim> pkfhat(nRefined,scaledReciprocalLattice);
            FFT::fft(pkf.values.template cast<dcomplex>(),pkfhat.values);

            // Correct pkfkhat as the kernel function was centered w.r.t the scaled lattice
            LatticeFunction<dcomplex,dim> temp(nRefined,scaledReciprocalLattice);
            /*
            Eigen::Vector<double,dim> center= scaledLattice.latticeBasis.col(0)/2 +
                                              scaledLattice.latticeBasis.col(1)/2 +
                                              scaledLattice.latticeBasis.col(2)/2;
                                              */
           Eigen::Vector<double,dim> center= scaledLattice.latticeBasis.rowwise().sum()/2;

            Exponential<dim> exp_function(center);
            LatticeFunction<dcomplex,dim> exp_factor(nRefined,scaledReciprocalLattice,exp_function);
            temp.values= exp_factor.values * pkfhat.values;
            pkfhat.values= temp.values;

            // Condense pkfhat to the unscaled lattice
            Lattice<dim> reciprocalLattice(m.reciprocalBasis);
            LatticeFunction<dcomplex,dim> pkfhatReduced(n,reciprocalLattice);
            Eigen::array<Eigen::DenseIndex, dim> strides;
            strides.fill(1);
            for(int i=0; i<dim; i++)
                strides[i] = factor(i);
            pkfhatReduced.values= pkfhat.values.stride(strides);

            // Fourier transform of the multilattice
            LatticeFunction<dcomplex,dim> mhat(n,reciprocalLattice);
            for(const auto& atom :  m.basisAtoms.colwise()) {
                Exponential<dim> exp_m(-atom);
                LatticeFunction<dcomplex,dim> temp(n,reciprocalLattice,exp_m);
                mhat.values+= temp.values;
            }

            FFT::ifft(mhat.values*pkfhatReduced.values,tempOutput.values);
            //if (typeid(Scalar) == typeid(dcomplex))
            //    return tempOutput;
            output.values = tempOutput.values.real();
            return output;
        }
    };

    template<typename Scalar, int dim>
    basic_ostream<char>& operator<<(basic_ostream<char>& s, const PeriodicFunction<Scalar, dim>& fun)
    {
        s << "x coord, y coord, z coord, scalar" << std::endl;
        LatticeVector<dim> x(fun.lattice);
        auto n= fun.values.dimensions();
        for (int i = 0; i < n[0]; i++) {
            for (int j = 0; j < n[1]; j++) {
                for (int k = 0; k < n[2]; k++) {
                    Eigen::Vector<double,dim> x= i*fun.lattice.latticeBasis.col(0)/n[0] +
                                                 j*fun.lattice.latticeBasis.col(1)/n[1] +
                                                 k*fun.lattice.latticeBasis.col(2)/n[2];
                    const Eigen::IOFormat fmt(15, 0, ", ", "", "\t", "", "", "");
                    s << x.transpose().format(fmt) << ", " << std::setw(25) << std::setprecision(15) << fun.values(i,j,k) << std::endl;
                }
            }
        }
        return s;
    }

}

#endif //OILAB_PERIODICFUNCTION_H
