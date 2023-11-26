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
    const double A1= -84.6841/factor;
    const double A2= 120.6472/factor;
    const double A3= 142.2017/factor;
    const double A4= -73.3189/factor;


    const double a_planar_1= 1.00316;
    const double a_planar_2= 2.51913;
    const double a_planar_3= 2.44581;
    const double a_planar_4= 0.56255;
    const double a_perp_1= 0.27752;
    const double a_perp_2= 2.29511;
    const double a_perp_3= 2.925396;
    const double a_perp_4= 1.32379;
public:
    explicit Mayer(double _domainSize) : Function<Mayer<dim>,double,dim>(_domainSize){}
    double operator()(const Eigen::Vector<double,dim>& vec) const
    {
        double output= 0.0;
        double exponent= a_planar_1*pow(vec(0),2) + a_planar_1*pow(vec(1),2) + a_perp_1*pow(vec(2),2);
        output= A1*exp(-exponent);
        exponent= a_planar_2*pow(vec(0),2) + a_planar_2*pow(vec(1),2) + a_perp_2*pow(vec(2),2);
        output= output + A2*exp(-exponent);
        exponent= a_planar_3*pow(vec(0),2) + a_planar_3*pow(vec(1),2) + a_perp_3*pow(vec(2),2);
        output= output + A3*exp(-exponent);
        exponent= a_planar_4*pow(vec(0),2) + a_planar_4*pow(vec(1),2) + a_perp_4*pow(vec(2),2);
        output= output + A4*exp(-exponent);
        return output;
    }
};
#endif //OILAB_MAYER_H
