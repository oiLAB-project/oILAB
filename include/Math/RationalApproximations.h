//
// Created by Nikhil Chandra Admal on 4/20/23.
//

#ifndef OILAB_RATIONALAPPROXIMATIONS_H
#define OILAB_RATIONALAPPROXIMATIONS_H

#include <IntegerMath.h>
#include <Rational.h>
#include <Farey.h>

namespace gbLAB {
    template<typename IntScalarType>
    class RationalApproximations
    {
    public:
        double number;
        double tolerance;
        IntScalarType max_denominator;
        std::vector<Rational<IntScalarType>> approximations;

        RationalApproximations(double number, int max_denominator, double tolerance):
        /* init */ number(number),
        /* init */ max_denominator(max_denominator),
        /* init */ tolerance(tolerance)
        {
            double n=floor(number);
            auto seq = farey(max_denominator,false);
            for (auto p : seq) {
                double approx = static_cast<double>(p.first) / static_cast<double>(p.second);
                if (std::abs(number - n - approx) <= tolerance) {
                    p.first= p.first + n*p.second;
                    approximations.push_back(Rational<IntScalarType>(p.first,p.second));
                }
            }
        }




    };

}
#endif //OILAB_RATIONALAPPROXIMATIONS_H
