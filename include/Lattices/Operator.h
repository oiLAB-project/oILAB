//
// Created by Nikhil Chandra Admal on 7/5/23.
//

#ifndef OILAB_OPERATOR_H
#define OILAB_OPERATOR_H

#include <LatticeModule.h>
#include <unsupported/Eigen/CXX11/Tensor>

namespace gbLAB {
    template<typename Derived, int dim>
    class Operator {
    public:
        using Scalar = double;
        const Derived &derivedOperator;
        const Lattice<dim> L;
        const Eigen::array<Eigen::Index, dim> n;

        Operator(const Eigen::Matrix<double,dim,dim>& A,
                 const Eigen::array<Eigen::Index,dim>& n_) : derivedOperator(static_cast<const Derived &>(*this)),
                                                             L(A),
                                                             n(n_)
        {
        }

        auto domain() const { return L.latticeBasis; }

        Eigen::Index rows() const { return std::accumulate(begin(n), end(n), 1, std::multiplies<>()); }
        Eigen::Index cols() const { return std::accumulate(begin(n), end(n), 1, std::multiplies<>()); }

        void perform_op(const double *x_in, double *y_out) const {
            derivedOperator.perform_op(x_in, y_out);
        }

        // implement expression templates for operators

    };

    template<typename E1, typename E2, int dim>
    class OperatorSum : public Operator<OperatorSum<E1,E2,dim>,dim> {
        const E1& o1;
        const E2& o2;
    public:
        OperatorSum(const E1 &o1_, const E2 &o2_) : Operator<OperatorSum<E1,E2,dim>,dim>(o1_.L.latticeBasis,o1_.n),
                                                    o1(o1_), o2(o2_)

        {
            // ensure that the domains of the two operators and their discretizations match
            assert(o1.n == o2.n && o1.domain().isApprox(o2.domain()));
        }

        void perform_op(const double *x_in, double *y_out) const {
            int nx = this->rows();
            Eigen::VectorXd y1(nx);
            Eigen::VectorXd y2(nx);
            o1.perform_op(x_in, y1.data());
            o2.perform_op(x_in, y2.data());
            Eigen::Map<Eigen::VectorXd> y(y_out, nx);
            y = y1 + y2;
        }
    };

    template<typename E1, typename E2, int dim>
    OperatorSum<E1, E2, dim> operator+(const Operator<E1, dim>& u, const Operator<E2,dim>& v) {
        return OperatorSum<E1, E2, dim>(*static_cast<const E1 *>(&u), *static_cast<const E2 *>(&v));
    }

    template<typename T, typename E, int dim>
    class OperatorScalarProduct : public Operator<OperatorScalarProduct<T,E,dim>,dim> {
        double s;
        const E &op;
    public:
        OperatorScalarProduct(const T& s_, const E& op_) : Operator<OperatorScalarProduct<T,E,dim>,dim>(op_.L.latticeBasis,op_.n),
                                                           s(s_), op(op_)

        {}

        void perform_op(const double* x_in, double* y_out) const {
            int nx = this->rows();
            op.perform_op(x_in, y_out);
            Eigen::Map<Eigen::VectorXd> y(y_out, nx);
            y = s * y;
        }
    };

    template<typename T, typename E, int dim>
    OperatorScalarProduct<T, E, dim> operator*(const T &s, const Operator<E, dim> &u) {
        return OperatorScalarProduct<T, E, dim>(s, *static_cast<const E *>(&u));
    }

}
#endif //OILAB_OPERATOR_H
