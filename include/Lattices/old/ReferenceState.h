//
// Created by Nikhil Chandra Admal on 3/27/24.
//

#ifndef OILAB_REFERENCESTATE_H
#define OILAB_REFERENCESTATE_H

#include "Gb.h"
#include "OrderedTuplet.h"

namespace gbLAB {

    template<int dim>
    class ReferenceState
    {
        using IntScalarType= typename LatticeCore<dim>::IntScalarType;
    public:
        const Gb<dim>& gb;
        const ReciprocalLatticeVector<dim>& axis;
        const int periodScaling;

        std::map<Triplet> refState;
        Eigen::VectorXd planeEnergies;
        //std::map<LatticeVector<dim>,LatticeVector<dim>> coincidence;

        // optionally include the cell U
        explicit ReferenceState(const Gb<dim>& gb,
                                const ReciprocalLatticeVector<dim>& axis,
                                const int& periodScaling);
        int numberOfPlanesOrthogonalToGB() const;
        LatticeVector<dim> shiftLatticeVector() const;
    };



} // gbLAB

#endif //OILAB_REFERENCESTATE_H
