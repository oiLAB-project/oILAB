//
// Created by Nikhil Chandra Admal on 7/9/23.
//

#ifndef OILAB_POLYNOMIAL_H
#define OILAB_POLYNOMIAL_H

#include <deque>
#include <cmath>

template<typename Scalar>
class Polynomial {
public:
    Polynomial(const std::deque<Scalar>& coefficients):coefficients(coefficients)
    {}
    Polynomial()=default;
    virtual ~Polynomial()= default;
    std::deque<Scalar> coefficients;

    Scalar operator()(const double& argument) const
    {
        int degree= 0;
        Scalar value= 0;
        for(const auto& coefficient : coefficients)
        {
            value= value + coefficient*pow(argument,degree);
            degree++;
        }

        return value;
    }

};

template<typename Derived, typename Scalar, int dim>
class Function
{
private:
    const Derived& derivedFunction;
public:
    Function():derivedFunction(static_cast<const Derived&>(*this))
    {}
    Scalar operator()(const Eigen::Vector<double,dim>& vec) const
    {
        return derivedFunction(vec);
    }
};

template<typename Scalar, int dim>
class PiecewisePolynomial : public Function<PiecewisePolynomial<Scalar,dim>,Scalar,dim>
{
private:
    std::map<std::pair<double,double>,Polynomial<Scalar>> piecewisePolynomial;

public:
    explicit PiecewisePolynomial(const std::map<double,Scalar>& map)
    {
        // ensure the domain of map is non-negative
        assert(  (map.cbegin())->first >= 0  );

        for (auto iter = map.cbegin(); iter != map.cend(); iter++)
        {
            // (y-y1)/(y2-y1) = (x-x1)/(x2-x1)
            // y = (y2-y1)/(x2-x1)  *  (x-x1)   +   y1
            //   = y1-x1*(y2-y1)/(x2-x1) + (y2-y1)/(x2-x1)  *   x
            auto nextiter= std::next(iter);
            if (nextiter != map.cend())
            {
                double x1= iter->first;  double x2= nextiter->first;
                std::pair<double,double> interval{x1,x2};

                Scalar y1= iter->second; Scalar y2= nextiter->second;
                std::deque<Scalar> linear{ y1-x1*(y2-y1)/(x2-x1),  (y2-y1)/(x2-x1) };
                Polynomial<Scalar> linearPolynomial{linear};
                piecewisePolynomial[interval]= linearPolynomial;
            }
        }
    }

    Scalar operator()(const Eigen::Vector<double,dim>& vec) const
    {
        double r= vec.norm();

        auto iter = std::find_if(piecewisePolynomial.cbegin(), piecewisePolynomial.cend(),
                                 [=](const std::pair< std::pair<double,double>, Polynomial<Scalar>>& polynomial)
                                 {
                                     return r>=polynomial.first.first &&  r<polynomial.first.second;
                                 });

        if (iter == piecewisePolynomial.end()) return 0.0;
        else
            return (iter->second)(r);
    }
};

template<int dim>
class Exponential : public Function<Exponential<dim>,std::complex<double>,dim>
{
private:
    Eigen::Vector<double,dim> x;
public:
    Exponential(Eigen::Vector<double,dim> _x) : x(_x)
    {}
    std::complex<double> operator()(const Eigen::Vector<double,dim>& vec) const
    {
        return exp(2*M_PI*std::complex<double>(0,1)*vec.dot(x));
    }
};

#endif //OILAB_POLYNOMIAL_H
