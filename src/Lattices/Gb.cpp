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
    /* init */,nA(&n.lattice == &(bc.A) ? bc.A.reciprocalLatticeDirection(n.cartesian()) : bc.A.reciprocalLatticeDirection(-n.cartesian()))
    /* init */,nB(&n.lattice == &(bc.B) ? bc.B.reciprocalLatticeDirection(n.cartesian()) : bc.B.reciprocalLatticeDirection(-n.cartesian()))
    {
    }

    template<int dim>
    double Gb<dim>::stepHeightA(const LatticeVector<dim>& d) const
    {
        ReciprocalLatticeDirection<dim> dir= bc.getReciprocalLatticeDirectionInC(nA.reciprocalLatticeVector());
        double cslPlaneSpacing= dir.planeSpacing();
        double step= bc.shiftTensorA(d).cartesian().dot(nA.cartesian().normalized());
        return std::remainder(step,cslPlaneSpacing);
    }

    template<int dim>
    double Gb<dim>::stepHeightB(const LatticeVector<dim>& d) const
    {
        ReciprocalLatticeDirection<dim> dir= bc.getReciprocalLatticeDirectionInC(nA.reciprocalLatticeVector());
        double cslPlaneSpacing= dir.planeSpacing();
        double step= bc.shiftTensorB(d).cartesian().dot(nB.cartesian().normalized());
        return std::remainder(step,cslPlaneSpacing);
    }

    template class Gb<2>;
    template class Gb<3>;
    template class Gb<4>;
    template class Gb<5>;
}


#endif //OILAB_GB_CPP
