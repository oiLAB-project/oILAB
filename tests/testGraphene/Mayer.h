//
// Created by Nikhil Chandra Admal on 7/19/23.
//

#ifndef OILAB_MAYER_H
#define OILAB_MAYER_H

#include <Function.h>
template<int dim>
class Mayer: public Function<Mayer<dim>,double,dim>
{
private:
    const double factor= 27.2114079527;
    const double A1= 10.607/factor;
    const double A2= 29.711/factor;
    const double A3= -98.911/factor;
    const double a1= 0.12126;
    const double a2= 1.9148;
    const double a3= 0.60078;
public:
    explicit Mayer(double _domainSize) : Function<Mayer<dim>,double,dim>(_domainSize){}
    double operator()(const Eigen::Vector<double,dim>& vec) const
    {
        double r2= vec.squaredNorm();
        return A1*exp(-a1*r2) + A2*exp(-a2*r2) + A3*exp(-a3*r2);
    }
};
#endif //OILAB_MAYER_H
