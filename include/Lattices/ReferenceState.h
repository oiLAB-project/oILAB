//
// Created by Nikhil Chandra Admal on 3/27/24.
//

#ifndef OILAB_REFERENCESTATE_H
#define OILAB_REFERENCESTATE_H

#include <Gb.h>
#include <Triplet.h>

namespace gbLAB {

    template<int dim>
    class ReferenceState
    {
        using IntScalarType= typename LatticeCore<dim>::IntScalarType;
    public:
        const Gb<dim>& gb;
        const ReciprocalLatticeVector<dim>& axis;
        const int periodScaling;
        std::vector<Triplet> refState;
        Eigen::VectorXd planeEnergies;

        explicit ReferenceState(const Gb<dim>& gb,
                                const ReciprocalLatticeVector<dim>& axis,
                                const int& periodScaling);
        int numberOfPlanesOrthogonalToGB() const;
        LatticeVector<dim> shiftLatticeVector() const;
    };



} // gbLAB

#endif //OILAB_REFERENCESTATE_H
