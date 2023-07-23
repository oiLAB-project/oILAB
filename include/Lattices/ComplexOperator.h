//
// Created by Nikhil Chandra Admal on 7/15/23.
//

#ifndef OILAB_DIFFCOMPLEX_H
#define OILAB_DIFFCOMPLEX_H


#include <Operator.h>
#include <LatticeModule.h>
#include <unsupported/Eigen/CXX11/Tensor>
#include <FFT.h>


namespace gbLAB {
    template<typename T1, typename T2, int dim>
    class ComplexOperator : public Operator<ComplexOperator<T1,T2,dim>,dim>
    {
    public:
        using Scalar= double;
        using Derived=  Operator<ComplexOperator<T1,T2,dim>,dim>;
        const Operator<T1,dim>& opReal;
        const Operator<T2,dim>& opImag;

        explicit ComplexOperator(const T1& opReal_, const T2& opImag_) :
                        Derived(opReal_.L.latticeBasis, opReal_.n),
                               /*
                               [&opReal_]()->Eigen::array<Eigen::Index,dim>{
                                    Eigen::array<Eigen::Index,dim> temp;
                                    for(int i=0; i<dim; i++)
                                        temp[i]=2*opReal_.n[i];
                                    return temp;
                                }()),
                                */
                               opReal(opReal_), opImag(opImag_)
        {
            assert(opReal.rows() == opImag.rows() && opReal.cols() == opImag.cols() && opReal.domain().isApprox(opImag.domain()));
        }

        Eigen::Index rows() const { return 2*std::accumulate(begin(this->n), end(this->n), 1, std::multiplies<>()); }
        Eigen::Index cols() const { return 2*std::accumulate(begin(this->n), end(this->n), 1, std::multiplies<>()); }

        void perform_op(const double* x_in, double* y_out) const
        {
            /* A = B + i C
             *
             * D = [B   -C]
             *     [C    B]
             * A and D have the same eigenvalues
             *
             */
            int nx = rows();
            Eigen::Map<const Eigen::VectorXd> x_in_topHalf(x_in,nx/2);
            Eigen::Map<const Eigen::VectorXd> x_in_bottomHalf(x_in+nx/2,nx/2);
            Eigen::Map<Eigen::VectorXd> y_out_topHalf(y_out,nx/2);
            Eigen::Map<Eigen::VectorXd> y_out_bottomHalf(y_out+nx/2,nx/2);

            Eigen::VectorXd temp(nx/2);
            // B xt
            opReal.perform_op(x_in_topHalf.data(),y_out_topHalf.data());
            // C xb
            opImag.perform_op(x_in_bottomHalf.data(),temp.data());
            // B*xt - C*xb
            y_out_topHalf= y_out_topHalf - temp;


            // C xt
            opImag.perform_op(x_in_topHalf.data(),y_out_bottomHalf.data());
            // B xb
            opReal.perform_op(x_in_bottomHalf.data(),temp.data());
            // C*xt + B*xb
            y_out_bottomHalf= y_out_bottomHalf+ temp;
        }
    };
}

#endif //OILAB_DIFFCOMPLEX_H
