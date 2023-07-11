//
// Created by Nikhil Chandra Admal on 7/5/23.
//

#ifndef OILAB_OPERATOR_H
#define OILAB_OPERATOR_H

#include <Spectra/GenEigsSolver.h>

template<typename Derived>
class Operator
{
public:
    using Scalar = double;
    const Derived& derivedOperator;

    Operator():derivedOperator(static_cast<const Derived&>(*this))
    {}

    auto domain() const
    {
        return derivedOperator.L.latticeBasis;
    }
    Eigen::Index rows() const
    {
        return derivedOperator.rows();
    }
    Eigen::Index cols() const
    {
        return derivedOperator.cols();
    }
    void perform_op(const double* x_in, double* y_out) const
    {
        derivedOperator.perform_op(x_in,y_out);
    }

    // implement expression templates for operators

};

template<typename E1, typename E2>
class OperatorSum : public Operator<OperatorSum<E1,E2>>
{
    const E1& o1;
    const E2& o2;
public:
    OperatorSum(const E1& o1_, const E2& o2_) : o1(o1_),o2(o2_)
    {
        // ensure that the domains of the two operators and their discretizations match
        assert(o1.rows() == o2.rows() && o1.cols() == o2.cols() && o1.domain().isApprox(o2.domain()));
    }
    OperatorSum(const E1&& o1_, const E2&& o2_) : o1(o1_),o2(o2_)
    {
        // ensure that the domains of the two operators and their discretizations match
        assert(o1.rows() == o2.rows() && o1.cols() == o2.cols() && o1.domain().isApprox(o2.domain()));
    }
    Eigen::Index rows() const
    {
        return o1.rows();
    }
    Eigen::Index cols() const
    {
        return o1.cols();
    }
    void perform_op(const double* x_in, double* y_out) const
    {
        int nx = rows();
        Eigen::VectorXd y1(nx);
        Eigen::VectorXd y2(nx);
        o1.perform_op(x_in,y1.data());
        o2.perform_op(x_in,y2.data());
        Eigen::Map<Eigen::VectorXd> y(y_out,nx);
        y= y1+y2;
    }
};

template <typename E1, typename E2>
OperatorSum<E1, E2> operator+(const Operator<E1>& u, const Operator<E2>& v)
{
    return OperatorSum<E1, E2>(*static_cast<const E1*>(&u), *static_cast<const E2*>(&v));
}

template <typename E1, typename E2>
OperatorSum<E1, E2> operator+(const Operator<E1>&& u, const Operator<E2>&& v)
{
    return OperatorSum<E1, E2>(static_cast<const E1&&>(*static_cast<const E1*>(&u)),
                               static_cast<const E2&&>(*static_cast<const E2*>(&v)));
}



#endif //OILAB_OPERATOR_H
