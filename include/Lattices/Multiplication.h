//
// Created by Nikhil Chandra Admal on 7/6/23.
//

#ifndef OILAB_MULTIPLICATION_H
#define OILAB_MULTIPLICATION_H

#include <Operator.h>
#include <LatticeModule.h>
#include <PeriodicFunction.h>

namespace gbLAB {
    template<int dim>
    class Multiplication : public Operator<Multiplication<dim>,dim> {
    public:
        const PeriodicFunction<double,dim>& pf;

        explicit Multiplication(const PeriodicFunction<double,dim>& _pf) :
                Operator<Multiplication<dim>,dim>(_pf.lattice.latticeBasis,_pf.values.dimensions()),
                pf(_pf)
        { }

        void perform_op(const double* x_in, double* y_out) const
        {
            int n= this->rows();
            Eigen::Map<const Eigen::VectorXd> x(x_in,n);
            Eigen::Map<Eigen::VectorXd> y(y_out,n);
            Eigen::Map<const Eigen::VectorXd> data(pf.values.data(),n);
            y=x.cwiseProduct(data);
        }
    };
}
#endif //OILAB_MULTIPLICATION_H
