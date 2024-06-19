//
// Created by Nikhil Chandra Admal on 6/17/24.
//
#include <unsupported/Eigen/CXX11/Tensor>
#include <LatticeFunction.h>
#include <list>

#ifndef OILAB_MATERIALTENSORS_H
#define OILAB_MATERIALTENSORS_H

namespace gbLAB {
    class GbMaterialTensors {
    public:
        static double lambda, mu;
        static double tensorC(const int& k, const int& p, const int& l, const int& q);
        static std::complex<double> tensorFhat(const int& k, const int& l, const int& i, const int& j, const Eigen::Vector3d& xi);
        static std::complex<double> tensorGhat(const int& i, const int& k, const int& t, const int& r, const Eigen::VectorXd& xi);
        static std::complex<double> tensorHhat(const int& t, const int& i, const Eigen::VectorXd &xi);
    };

}


#endif //OILAB_MATERIALTENSORS_H
