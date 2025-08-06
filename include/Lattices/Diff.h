//
// Created by Nikhil Chandra Admal on 7/14/23.
//

#ifndef OILAB_DIFF_H
#define OILAB_DIFF_H

#include <Operator.h>
#include <LatticeModule.h>
#include <unsupported/Eigen/CXX11/Tensor>
#include <FFT.h>


namespace gbLAB {
    template<int dim>
    class Diff : public Operator<Diff<dim>,dim>
    {
    public:
        using dcomplex= std::complex<double>;
        using Operator<Diff<dim>,dim>::n, Operator<Diff<dim>,dim>::L;
        const Eigen::array<Eigen::Index,dim> d;

        explicit Diff(const Eigen::array<Eigen::Index,dim>& d_,
                      const Eigen::Matrix<double,dim,dim>& A,
                      const Eigen::array<Eigen::Index,dim>& n_) :
                      Operator<Diff<dim>,dim>(A,n_),
                      d(d_)
        {
            for(const auto& order : d)
                assert(order >= 0);

        }

        template<int dm=dim>
        typename std::enable_if<dm==1,void>::type
        perform_op(const double* x_in, double* y_out) const
        {
            Eigen::TensorMap<const Eigen::Tensor<double,dm>> xReal(x_in,n);
            const Eigen::Tensor<dcomplex,dm> x=xReal.template cast<dcomplex>();
            Eigen::TensorMap<Eigen::Tensor<double,dm>> y(y_out,n);

            // if d=0 y_out= x_in and return
            int totalOrder= 0;
            for(const auto& order : d)
            {
                totalOrder= totalOrder + order;
            }
            if (totalOrder == 0) {
                y= xReal;
                return;
            }

            // Compute y = Lx using FFT
            Eigen::Tensor<dcomplex,dm> xhat(n);
            xhat.setZero();
            FFT::fft(x,xhat);

            Eigen::Tensor<dcomplex,dm> d2fhat(n);
            d2fhat.setZero();

            // only part which is dimension dependent
            for (int i = 0; i < n[0]; ++i) {
                ReciprocalLatticeVector<dm> r(L);
                dcomplex factor(1,0);
                for (int k= 0; k<dim; ++k) {
                    if (d[k] == 0) continue;
                    else if (d[k] % 2 == 0)
                        r << (i <= n[0] / 2 ? i : i - n[0]);
                    else
                        r << (i == n[0] / 2 ? 0 : -n[0] * (i / (n[0] / 2)) + i);
                    factor = factor * std::pow(-2.0 * M_PI *
                                               dcomplex(0, 1) *
                                               r.cartesian()(k),
                                               d[k]);
                }
                d2fhat(i)= xhat(i) * factor;
            }

            // Laplacian Lf
            Eigen::Tensor<dcomplex,dm> Lf(n);
            Lf.setZero();
            FFT::ifft(d2fhat,Lf);
            y= Lf.real();
        }

        template<int dm=dim>
        typename std::enable_if<dm==2,void>::type
        perform_op(const double* x_in, double* y_out) const
        {
            Eigen::TensorMap<const Eigen::Tensor<double,2>> xReal(x_in,n);
            const Eigen::Tensor<dcomplex,2> x=xReal.cast<dcomplex>();
            Eigen::TensorMap<Eigen::Tensor<double,2>> y(y_out,n);

            // if d=0 y_out= x_in and return
            int totalOrder= 0;
            for(const auto& order : d)
            {
                totalOrder= totalOrder + order;
            }
            if (totalOrder == 0) {
                y= xReal;
                return;
            }

            // Compute y = Lx using FFT
            Eigen::Tensor<dcomplex,2> xhat(n);
            xhat.setZero();
            FFT::fft(x,xhat);

            Eigen::Tensor<dcomplex,2> d2fhat(n);
            d2fhat.setZero();

            // only part which is dimension dependent
                for (int i = 0; i < n[0]; ++i) {
                    for (int j = 0; j < n[1]; ++j) {
                        ReciprocalLatticeVector<2> r(L);
                        dcomplex factor(1,0);
                        for (int k= 0; k<dim; ++k) {
                            if (d[k] == 0) continue;
                            else if (d[k] % 2 == 0)
                                r << (i <= n[0] / 2 ? i : i - n[0]),
                                     (j <= n[1] / 2 ? j : j - n[1]);
                            else
                                r << (i == n[0] / 2 ? 0 : -n[0] * (i / (n[0] / 2)) + i),
                                     (j == n[1] / 2 ? 0 : -n[1] * (j / (n[1] / 2)) + j);
                            factor = factor * std::pow(-2.0 * M_PI *
                                                       dcomplex(0, 1) *
                                                       r.cartesian()(k),
                                                       d[k]);
                        }
                        d2fhat(i,j)= xhat(i,j) * factor;
                    }
                }

            // Laplacian Lf
            Eigen::Tensor<dcomplex,2> Lf(n);
            Lf.setZero();
            FFT::ifft(d2fhat,Lf);
            y= Lf.real();
        }

        template<int dm=dim>
        typename std::enable_if<dm==3,void>::type
        perform_op(const double* x_in, double* y_out) const
        {
            Eigen::TensorMap<const Eigen::Tensor<double,dm>> xReal(x_in,n);
            const Eigen::Tensor<dcomplex,dm> x=xReal.template cast<dcomplex>();
            Eigen::TensorMap<Eigen::Tensor<double,dm>> y(y_out,n);

            // if d=0 y_out= x_in and return
            int totalOrder= 0;
            for(const auto& order : d)
            {
                totalOrder= totalOrder + order;
            }
            if (totalOrder == 0) {
                y= xReal;
                return;
            }

            // Compute y = Lx using FFT
            Eigen::Tensor<dcomplex,dm> xhat(n);
            xhat.setZero();
            FFT::fft(x,xhat);

            Eigen::Tensor<dcomplex,dm> d2fhat(n);
            d2fhat.setZero();

            // only part which is dimension dependent
            for (int i = 0; i < n[0]; ++i) {
                for (int j = 0; j < n[1]; ++j) {
                    for (int k = 0; k < n[2]; ++k) {
                        ReciprocalLatticeVector<dm> r(L);
                        dcomplex factor(1, 0);
                        for (int l = 0; l < dm; ++l) {
                            if (d[l] == 0) continue;
                            else if (d[l] % 2 == 0)
                                r << (i <= n[0] / 2 ? i : i - n[0]),
                                     (j <= n[1] / 2 ? j : j - n[1]),
                                     (k <= n[2] / 2 ? k : k - n[2]);
                            else
                                r << (i == n[0] / 2 ? 0 : -n[0] * (i / (n[0] / 2)) + i),
                                     (j == n[1] / 2 ? 0 : -n[1] * (j / (n[1] / 2)) + j),
                                     (k == n[2] / 2 ? 0 : -n[2] * (k / (n[2] / 2)) + k);
                            factor = factor * std::pow(-2.0 * M_PI *
                                                       dcomplex(0, 1) *
                                                       r.cartesian()(l),
                                                       d[l]);
                        }
                        d2fhat(i,j,k) = xhat(i,j,k) * factor;
                    }
                }
            }

            // Laplacian Lf
            Eigen::Tensor<dcomplex,dm> Lf(n);
            Lf.setZero();
            FFT::ifft(d2fhat,Lf);
            y= Lf.real();
        }
    };
}

#endif //OILAB_DIFF_H
