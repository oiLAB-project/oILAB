//
// Created by Nikhil Chandra Admal on 7/6/23.
//

#ifndef OILAB_MULTIPLICATION_H
#define OILAB_MULTIPLICATION_H

#include <Operator.h>
#include <LatticeModule.h>
#include <PeriodicFunction.h>

namespace gbLAB {
    template<typename Scalar, int dim>
    class Multiplication : public Operator<Multiplication<Scalar,dim>> {
    private:
        int n;

    public:
        Lattice<dim> L;
        const PeriodicFunction<Scalar,dim>& pf;

        explicit Multiplication(const PeriodicFunction<Scalar,dim>& _pf) : pf(_pf), L(_pf.lattice)
        {
            n= std::accumulate(begin(pf.values.dimensions()), end(pf.values.dimensions()), 1, std::multiplies<>());
        }


        Eigen::Index rows() const { return n;}
        Eigen::Index cols() const { return n;}

        void perform_op(const double* x_in, double* y_out) const
        {
            Eigen::Map<const Eigen::VectorXd> x(x_in,n);
            Eigen::Map<Eigen::VectorXd> y(y_out,n);
            Eigen::Map<const Eigen::VectorXd> data(pf.values.data(),n);
            y=x.cwiseProduct(data);
        }
    };
}
#endif //OILAB_MULTIPLICATION_H
