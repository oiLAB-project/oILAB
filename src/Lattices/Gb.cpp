//
// Created by Nikhil Chandra Admal on 11/5/22.
//
#ifndef OILAB_GB_CPP
#define OILAB_GB_CPP

#include "Gb.h"

namespace gbLAB
{
    template<int dim>
    Gb<dim>::Gb(const BiCrystal<dim>& bc, const ReciprocalLatticeDirection<dim>& n) :
    /* init */ bc(bc)
    /* init */,nA(&n.lattice == &(bc.A) ? bc.getReciprocalLatticeDirectionInA(n.reciprocalLatticeVector()) : bc.getReciprocalLatticeDirectionInA(-1*n.reciprocalLatticeVector()))
    /* init */,nB(&n.lattice == &(bc.B) ? bc.getReciprocalLatticeDirectionInB(n.reciprocalLatticeVector()) : bc.getReciprocalLatticeDirectionInB(-1*n.reciprocalLatticeVector()))
    {
    }

    template<int dim>
    double Gb<dim>::stepHeightA(const LatticeVector<dim>& d) const
    {
        return bc.shiftTensorA(d).cartesian().dot(nA.cartesian().normalized());
    }

    template<int dim>
    double Gb<dim>::stepHeightB(const LatticeVector<dim>& d) const
    {
        return bc.shiftTensorB(d).cartesian().dot(nB.cartesian().normalized());
    }

    template class Gb<2>;
    template class Gb<3>;
    template class Gb<4>;
    template class Gb<5>;
}


#endif //OILAB_GB_CPP
