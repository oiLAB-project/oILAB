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
#include <OrderedTuplet.h>

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

        static GbLatticeFunctions HhatInvComponents;
        //static FunctionFFTPair pipihat;
        static std::map<OrderedTuplet<dim>,PeriodicFunction<double, dim - 1>> piPeriodicFunctions;
        static std::map<OrderedTuplet<dim>,LatticeFunction<std::complex<double>, dim - 1>> pihatLatticeFunctions;
        // GBMesostateEnsemble should generate the bicrystal (member variable <OrderedTuplet,VectorDimD>) and pass it as a reference to each mesostate
        // pipihat should be map from OrderedTuplet to FunctionFFTPair. should be computed once in calculateb
        // at the same time, compute pipihat once
        // change xuPairs type to <Tiplet,VectorDimD>

        FunctionFFTPair bbhat;
        static FunctionFFTPair calculateb(const Eigen::Matrix<double, dim,dim-1>& domain,
                                          const std::map<OrderedTuplet<dim>,VectorDimD>& xuPairs,
                                          const std::array<Eigen::Index,dim-1>& n,
                                          const std::map<OrderedTuplet<dim>,VectorDimD>& points);
        static GbLatticeFunctions getHhatInvComponents(const Eigen::Matrix<double, dim,dim-1>& domain,
                                                       const std::array<Eigen::Index,dim-1>& n);
        static PeriodicFunction<double,dim-1>get_pi(const Eigen::Matrix<double,dim,dim-1>& domain,
                                                    const std::array<Eigen::Index,dim-1>& n,
                                                    const VectorDimD& point);
        static LatticeFunction<std::complex<double>,dim-1>get_pihat(const Eigen::Matrix<double,dim,dim-1>& domain,
                                                                    const std::array<Eigen::Index,dim-1>& n,
                                                                    const VectorDimD& point);

    public:


        const Eigen::Matrix<double,dim,dim-1> gbDomain;
        const std::map<OrderedTuplet<dim>,VectorDimD> xuPairs;
        std::array<Eigen::Index,dim-1> n;
        std::vector<PeriodicFunction<double,dim-1>> b;
        std::vector<LatticeFunction<std::complex<double>,dim-1>> bhat;
        std::map<OrderedTuplet<dim>,VectorDimD> atoms;

        double energy;

        GbContinuum(const Eigen::Matrix<double, dim,dim-1>& domain,
                    const std::map<OrderedTuplet<dim>, VectorDimD>& xuPairs,
                    const std::array<Eigen::Index,dim-1>& n,
                    const std::map<OrderedTuplet<dim>,VectorDimD>& atoms);

        VectorDimD displacement(const OrderedTuplet<dim>& x) const;
    };


    template<int dim>
    std::vector<LatticeFunction<std::complex<double>,dim-1>> GbContinuum<dim>::HhatInvComponents;

    template<int dim>
    std::map<OrderedTuplet<dim>,PeriodicFunction<double, dim - 1>> GbContinuum<dim>::piPeriodicFunctions;

    template<int dim>
    std::map<OrderedTuplet<dim>,LatticeFunction<std::complex<double>, dim - 1>> GbContinuum<dim>::pihatLatticeFunctions;
}
#endif //OILAB_GBPLASTICITY_H
