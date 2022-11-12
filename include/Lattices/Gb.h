//
// Created by Nikhil Chandra Admal on 11/5/22.
//

#ifndef OILAB_GB_H
#define OILAB_GB_H

#include "BiCrystal.h"
#include "ReciprocalLatticeDirection.h"

namespace gbLAB
{
    template<int dim>
    class Gb
    {
    using VectorDimI = typename LatticeCore<dim>::VectorDimI ;

    public:
        const BiCrystal<dim>& bc;
        const ReciprocalLatticeDirection<dim> nA;
        const ReciprocalLatticeDirection<dim> nB;
//        const double cslPlaneSpacing;
//        const double dsclPlaneSpacing;

        double stepHeightA(const LatticeVector<dim>& d) const;
        double stepHeightB(const LatticeVector<dim>& d) const;

        Gb(const BiCrystal<dim>& bc, const ReciprocalLatticeDirection<dim>& n);

    };
}

#endif //OILAB_GB_H
