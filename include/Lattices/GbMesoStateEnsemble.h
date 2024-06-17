//
// Created by Nikhil Chandra Admal on 5/17/24.
//

#ifndef OILAB_GBMESOSTATES_H
#define OILAB_GBMESOSTATES_H

#include <GbShifts.h>
#include <deque>
#include <GbMesoState.h>

namespace gbLAB {
    template<int dim>
    class GbMesoStateEnsemble : public GbShifts<dim>
    {
        using VectorDimD = typename LatticeCore<dim>::VectorDimD;
        static std::deque<GbMesoState<dim>> enumerateMesoStates(const Gb<dim>& gb,
                                                                const ReciprocalLatticeVector<dim>& axis,
                                                                const std::vector<std::pair<LatticeVector<dim>,VectorDimD>>& bShiftPairs);
    public:
        GbMesoStateEnsemble(const Gb<dim>& gb, const ReciprocalLatticeVector<dim>& axis);
        std::deque<GbMesoState<dim>> mesoStates;
    };

}

#endif //OILAB_GBMESOSTATES_H