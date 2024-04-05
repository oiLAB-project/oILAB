//
// Created by Nikhil Chandra Admal on 12/11/23.
//

#ifndef OILAB_MESOSTATE_H
#define OILAB_MESOSTATE_H

#include <ReferenceState.h>
#include <Gb.h>
#include <set>
#include <Dislocations.h>
#include <Triplet.h>

namespace gbLAB {

    template<int dim>
    class MesoState : public ReferenceState<dim>, public Dislocations
    {
        using IntScalarType= typename LatticeCore<dim>::IntScalarType;
        using VectorDimD = LatticeCore<3>::VectorDimD;
        using VectorDimI= typename LatticeCore<dim>::VectorDimI;
        using Matrix= typename Eigen::Matrix<IntScalarType,Eigen::Dynamic,Eigen::Dynamic>;
        using Vector2d= Eigen::Vector2d;
        using Matrix2d= Eigen::Matrix2d;

    public:
        // currently storing states of only one lattice as we are restricted to STGB
        std::vector<Triplet> currentState;
        std::vector<Triplet> defectsIndices;

        explicit MesoState(const Gb<dim>& gb,
                           const ReciprocalLatticeVector<dim>& axis,
                           const int& periodScaling,
                           const double& a2,
                           const int& nImages);

        explicit MesoState(const ReferenceState<dim>& rS, const double& a2, const int& nImages);

        double energy() const;
        Eigen::VectorXi getLocalStateCount(const int& numberOfInteractingPlanes) const
        {
            Eigen::VectorXi output(numberOfInteractingPlanes);
            for(int i=0; i<numberOfInteractingPlanes; ++i)
                output(i) = getOrthogonalPlaneIndices(i).size();
            return output;
        }

        void insertDislocation(const Triplet&);
        Triplet removeRandomDislocation();
        void removeDislocation(const Triplet& t);
        Triplet insertRandomDislocation(const int& dipoleSign);
        Triplet insertRandomDislocation();
        std::set<int> getOrthogonalPlaneIndices(const int& parallelPlaneIndex) const;

        //template<int dm=dim>
        typename std::enable_if<dim==3,void>::type
        box(const int& heightFactor,
            const int& dsclFactor,
            const std::string& name) const;

        template <typename T> int sgn(T val) const;
        bool operator<(const MesoState& rhs) const;
    };
}
#endif //OILAB_MESOSTATE_H
