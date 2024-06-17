//
// Created by Nikhil Chandra Admal on 5/26/24.
//

#ifndef OILAB_GBMESOSTATE_H
#define OILAB_GBMESOSTATE_H

#include <Gb.h>
#include <GbContinuum.h>
#include <LatticeCore.h>
namespace gbLAB {

    template<int dim>
    class GbMesoState : public GbContinuum<dim> {
        using VectorDimD = typename LatticeCore<dim>::VectorDimD;

        static Eigen::Matrix<double, dim,dim-1> getGbDomain(const std::vector<LatticeVector<dim>>& gbCslVectors);
        static std::deque<std::pair<VectorDimD, VectorDimD>> get_xuPairs(const std::deque<std::pair<LatticeVector<dim>, VectorDimD>>& bs,
                                                                         const VectorDimD& normal);
        static std::array<Eigen::Index,dim-1> discretize(const std::vector<LatticeVector<dim>>& gbCslVectors, const Gb<dim>& gb);

    public:

        const Gb<dim>& gb;
        const ReciprocalLatticeVector<dim>& axis;
        const std::vector<LatticeVector<dim>>& gbCslVectors;
        std::deque<std::pair<LatticeVector<dim>, VectorDimD>> bs;

        explicit GbMesoState(const Gb<dim> &gb,
                             const ReciprocalLatticeVector<dim>& axis,
                             const std::deque<std::pair<LatticeVector<dim>,VectorDimD>>& bs,
                             const std::vector<LatticeVector<dim>>& gbCslVectors);

        typename std::enable_if<dim==3,void>::type
        box(const int& heightFactor, const int& dsclFactor, const std::string& name) const;
    };
}
#endif //OILAB_GBMESOSTATE_H
