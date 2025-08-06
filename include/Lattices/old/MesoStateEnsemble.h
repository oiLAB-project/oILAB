//
// Created by Nikhil Chandra Admal on 3/25/24.
//

#ifndef OILAB_MESOSTATEENSEMBLE_H
#define OILAB_MESOSTATEENSEMBLE_H

#include"old/MesoState.h"
#include"old/ReferenceState.h"

namespace gbLAB {
    template <int dim>
    //class MesoStateEnsemble : public std::set<MesoState<dim>>
    class MesoStateEnsemble : public std::vector<MesoState<dim>>
    {
        using IntScalarType= typename LatticeCore<dim>::IntScalarType;

    public:
        //static int numberOfInteractingPlanes;
        const ReferenceState<dim>& rS;
        explicit MesoStateEnsemble(const ReferenceState<dim>& rS, const std::string& filename);
        explicit MesoStateEnsemble(const ReferenceState<dim>& rS, const int& numberOfMesoStates);
        explicit MesoStateEnsemble(const ReferenceState<dim>& rS);
        void read(const std::string& filename);
        void write(const std::string& filename) const;
        MesoStateEnsemble<dim> build() const;

    };
}

#endif //OILAB_MESOSTATEENSEMBLE_H
