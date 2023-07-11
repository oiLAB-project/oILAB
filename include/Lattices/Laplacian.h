//
// Created by Nikhil Chandra Admal on 7/5/23.
//

#ifndef OILAB_LAPLACIAN_H
#define OILAB_LAPLACIAN_H
#include <Operator.h>
#include <LatticeModule.h>
#include <unsupported/Eigen/CXX11/Tensor>
#include <FFT.h>


namespace gbLAB {
    template<int dim>
    class Laplacian : public Operator<Laplacian<dim>>
    {
    public:
        using dcomplex= std::complex<double>;
        Lattice<dim> L;
        Eigen::array<Eigen::Index,dim> n;

        explicit Laplacian(const Eigen::Matrix<double,dim,dim>& A, const Eigen::array<Eigen::Index,dim>& n_) : L(A), n(n_)
        { }
        Eigen::Index rows() const { return std::accumulate(begin(n), end(n), 1, std::multiplies<>()); }
        Eigen::Index cols() const { return std::accumulate(begin(n), end(n), 1, std::multiplies<>()); }

        template<int dm=dim>
        typename std::enable_if<dm==2,void>::type
        perform_op(const double* x_in, double* y_out) const
        {
            Eigen::TensorMap<const Eigen::Tensor<double,2>> xReal(x_in,n);
            const Eigen::Tensor<dcomplex,2> x=xReal.cast<dcomplex>();

            Eigen::TensorMap<Eigen::Tensor<double,2>> y(y_out,n);

            // Compute y = Lx using FFT
            Eigen::Tensor<dcomplex,2> xhat(n);
            xhat.setZero();
            FFT::fft(x,xhat);

            Eigen::Tensor<dcomplex,2> d2fhat(n);
            d2fhat.setZero();

            // only part which is dimension dependent
            for(int i=0; i<n[0]; ++i)
            {
                for (int j = 0; j < n[1]; ++j) {
                    ReciprocalLatticeVector<2> r(L);
                    r << (i <= n[0] / 2 ? i : i-n[0] ), (j <= n[1] / 2 ? j : j-n[1] );
                    d2fhat(i, j) = -xhat(i, j) * (dcomplex) (4.0 * std::pow(M_PI, 2) * r.cartesian().squaredNorm());
                }
            }

            // Laplacian Lf
            Eigen::Tensor<dcomplex,2> Lf(n);
            Lf.setZero();
            FFT::ifft(d2fhat,Lf);
            y= Lf.real();
        }
    };
}
#endif //OILAB_LAPLACIAN_H
