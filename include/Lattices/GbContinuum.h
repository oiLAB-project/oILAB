//
// Created by Nikhil Chandra Admal on 5/16/24.
//

#ifndef OILAB_GBPLASTICITY_H
#define OILAB_GBPLASTICITY_H

#include <LatticeCore.h>
#include "Eigen/Dense"
#include <PeriodicFunction.h>
#include <LatticeFunction.h>
#include <Function.h>

namespace gbLAB {

    class DisplacementKernel : public Function<DisplacementKernel, double> {
    private:
        Eigen::Vector<double, Eigen::Dynamic> normal;
    public:
        explicit DisplacementKernel(Eigen::Vector<double, Eigen::Dynamic> _normal, double domainSize) :
                normal(_normal), Function<DisplacementKernel, double>(domainSize)
        {}
        double operator()(const Eigen::Vector<double, Eigen::Dynamic> &x) const {
            return -x.dot(normal) / (4 * M_PI * (std::pow(x.norm(), 3)));
        }

    };

    class ShiftedDisplacementKernelFT :  public Function<ShiftedDisplacementKernelFT, std::complex<double>>{
    private:
        Eigen::Vector<double, Eigen::Dynamic> x, normal;
    public:
        explicit ShiftedDisplacementKernelFT(const Eigen::VectorXd& x,
                                             const Eigen::VectorXd& normal):
            x(x),normal(normal.normalized())
        { }
        std::complex<double> operator()(const Eigen::VectorXd& xi) const {
            double xNormalComponent=  x.dot(normal);
            auto xProjected= (x-xNormalComponent*normal).eval();
            double xiBar= xi.norm();
            std::complex<double> output= -0.5*std::exp(-2*M_PI* std::complex<double>(0,1) *
                                                       xProjected.dot(xi)
                                                     ) * std::exp(-2*M_PI*xiBar*abs(xNormalComponent));
            if (abs(xNormalComponent) < DBL_EPSILON)
                return output;
            else
                return output*xNormalComponent/abs(xNormalComponent);
        }

    };

    class H22InverseFT: public Function<H22InverseFT, std::complex<double>>
    {
    private:
        const Eigen::VectorXd e1,e3;
        const double lambda, mu;
    public:
        explicit H22InverseFT(const double& lambda,
                       const double& mu,
                       const Eigen::MatrixXd& domain) :
           lambda(lambda),
           mu(mu),
           e1(domain.col(0).normalized()),
           e3(domain.col(1).normalized())
        {}
        std::complex<double> operator()(const Eigen::VectorXd& xi) const
        {
            if (xi.isZero())
                return std::complex<double>(0,0);
            std::complex<double> output(0,0);
            double xi1= xi.dot(e1);
            double xi3= xi.dot(e3);
            double xiNorm= xi.norm();
            double nu= lambda/(2*(lambda+mu));
            double nuFactor= 1-nu;
            output= -(
                     // F1111
                    std::pow(lambda,2)* (
                        (std::pow(xi1,2)/std::pow(xiNorm,3)) -
                        3*std::pow(xi1,4)/(8*nuFactor*std::pow(xiNorm,5)) ) +
                    // F1212
                    (std::pow(lambda,2)-std::pow(mu,2))* (
                        -std::pow(xi1,2)/(8*nuFactor*std::pow(xiNorm,3)) ) +
                    // F1313
                    std::pow(lambda,2) * (
                        -std::pow(xi1,2)*std::pow(xi3,2)*3/(8*nuFactor*std::pow(xiNorm,5))
                            ) +
                    // F1122
                    (-lambda*mu-2*std::pow(mu,2)) * (
                            std::pow(xi1,2)/std::pow(xiNorm,3) -
                            std::pow(xi1,2)/(8*nuFactor*std::pow(xiNorm,3))
                            ) +
                    // F2211
                    std::pow(mu,2) * (
                            1./xiNorm -
                            std::pow(xi1,2)/(8*nuFactor*std::pow(xiNorm,3))
                            ) +
                    // F2222
                    (-lambda*mu-2*std::pow(mu,2)) * (
                            1./xiNorm -
                            3./(8*nuFactor*xiNorm)
                            ) +
                     // F2323
                     (-lambda*mu) * (
                             -std::pow(xi3,2)/(8*nuFactor*std::pow(xiNorm,3))
                             )
                     )/(8*M_PI*mu);
            return std::pow(output,-1);
        }
    };


    template<int dim>
    class GbContinuum {
        using VectorDimD = typename LatticeCore<dim>::VectorDimD;
        using FunctionFFTPair= typename std::pair<std::vector<PeriodicFunction<double,dim-1>>,
                                    std::vector<LatticeFunction<std::complex<double>,dim-1>>>;
    private:

        static VectorDimD H;
        static std::vector<LatticeFunction<std::complex<double>,dim-1>> pihats;
        static FunctionFFTPair calculateb(const Eigen::Matrix<double, dim,dim-1>& domain,
                                          const std::deque<std::pair<VectorDimD, VectorDimD>>& xuPairs,
                                          const std::array<Eigen::Index,dim-1>& n);
        static LatticeFunction<std::complex<double>,dim-1>
            calculateH22hatInverse(const Eigen::Matrix<double, dim,dim-1>& domain,
                            const std::array<Eigen::Index,dim-1>& n);
        static double calculateEnergy();
        FunctionFFTPair bbhat;


    public:

        static double lambda, mu;

        const Eigen::Matrix<double,dim,dim-1> gbDomain;
        const std::deque<std::pair<VectorDimD,VectorDimD>> xuPairs;
        std::array<Eigen::Index,dim-1> n;
        std::vector<PeriodicFunction<double,dim-1>> b;
        std::vector<LatticeFunction<std::complex<double>,dim-1>> bhat;
        double energy;

        GbContinuum(const Eigen::Matrix<double, dim,dim-1>& domain,
                    const std::deque<std::pair<VectorDimD, VectorDimD>>& xuPairs,
                    const std::array<Eigen::Index,dim-1>& n);

        VectorDimD displacement(const VectorDimD& x) const;
    };

template<int dim>
double GbContinuum<dim>::lambda;

template<int dim>
double GbContinuum<dim>::mu;

}
#endif //OILAB_GBPLASTICITY_H
