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
#include <GbMaterialTensors.h>

namespace gbLAB {

    class DisplacementKernel : public Function<DisplacementKernel, double> {
    private:
        Eigen::Vector<double, Eigen::Dynamic> normal;
    public:
        explicit DisplacementKernel(Eigen::Vector<double, Eigen::Dynamic> _normal, double domainSize);
        double operator()(const Eigen::Vector<double, Eigen::Dynamic> &x) const;
    };

    class ShiftedDisplacementKernelFT :  public Function<ShiftedDisplacementKernelFT, std::complex<double>>{
    private:
        Eigen::Vector<double, Eigen::Dynamic> x, normal;
    public:
        explicit ShiftedDisplacementKernelFT(const Eigen::VectorXd& x, const Eigen::VectorXd& normal);
        std::complex<double> operator()(const Eigen::VectorXd& xi) const;
    };

    class HhatInvFunction: public Function<HhatInvFunction, std::complex<double>>, GbMaterialTensors
    {
    private:
        const Eigen::VectorXd e1,e3;
    public:
        const int t,i;
        explicit HhatInvFunction(const int& t, const int& i, const Eigen::MatrixXd& domain);
        std::complex<double> operator()(const Eigen::VectorXd& xi) const;
    };

    template<int dim>
    class GbContinuum {
        using VectorDimD= typename LatticeCore<dim>::VectorDimD;
        using FunctionFFTPair= typename std::pair<std::vector<PeriodicFunction<double,dim-1>>,
                                                  std::vector<LatticeFunction<std::complex<double>,dim-1>>>;
        using GbLatticeFunctions= typename std::vector<LatticeFunction<std::complex<double>,dim-1>> ;
    private:

        static FunctionFFTPair calculateb(const Eigen::Matrix<double, dim,dim-1>& domain,
                                          const std::deque<std::pair<VectorDimD, VectorDimD>>& xuPairs,
                                          const std::array<Eigen::Index,dim-1>& n);

        static GbLatticeFunctions getHhatInvComponents(const Eigen::Matrix<double, dim,dim-1>& domain,
                                                       const std::array<Eigen::Index,dim-1>& n);

        static double calculateEnergy();
        static GbLatticeFunctions HhatInvComponents;
        FunctionFFTPair bbhat;

    public:


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
    std::vector<LatticeFunction<std::complex<double>,dim-1>> GbContinuum<dim>::HhatInvComponents;

}
#endif //OILAB_GBPLASTICITY_H
