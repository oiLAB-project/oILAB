//
// Created by Nikhil Chandra Admal on 7/19/23.
//

#ifndef OILAB_KP_H
#define OILAB_KP_H

#include <Function.h>
template<int dim>
class KP: public Function<KP<dim>,double,dim>
{
public:
    explicit KP(double _domainSize) : Function<KP<dim>,double,dim>(_domainSize){}
    double operator()(const Eigen::Vector<double,dim>& vec) const
    {
        double r2= vec.squaredNorm();
        if (r2<0.01)
            return 50.0;
        else
            return 0.0;
    }
};
#endif //OILAB_KP_H
