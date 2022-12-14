//
// Created by Nikhil Chandra Admal on 11/5/22.
//
#ifndef OILAB_GB_CPP
#define OILAB_GB_CPP

#include "Gb.h"

namespace gbLAB
{
    template<int dim>
    Gb<dim>::Gb(const BiCrystal<dim>& bc, const ReciprocalLatticeDirection<dim>& n)
    try :
    /* init */ bc(bc)
    /* init */,nA(&n.lattice == &(bc.A) ? n : bc.getReciprocalLatticeDirectionInA(-1*n.reciprocalLatticeVector()))
    /* init */,nB(&n.lattice == &(bc.B) ? n : bc.getReciprocalLatticeDirectionInB(-1*n.reciprocalLatticeVector()))
   // /* init */,nA(&n.lattice == &(bc.A) ? bc.A.reciprocalLatticeDirection(n.cartesian()) : bc.A.reciprocalLatticeDirection(-n.cartesian()))
   // /* init */,nB(&n.lattice == &(bc.B) ? bc.B.reciprocalLatticeDirection(n.cartesian()) : bc.B.reciprocalLatticeDirection(-n.cartesian()))
    {
        if (&n.lattice != &(bc.A) && &n.lattice != &(bc.B))
            throw std::runtime_error("The normal does not belong to the reciprocal lattices of A and B  \n");
    }
    catch(std::runtime_error& e)
    {
        std::cout << e.what() << std::endl;
        throw(std::runtime_error("GB construction failed. "));
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
