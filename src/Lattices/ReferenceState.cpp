//
// Created by Nikhil Chandra Admal on 3/27/24.
//

#include "ReferenceState.h"

namespace gbLAB {
    template<int dim>
    ReferenceState<dim>::ReferenceState(const Gb<dim>& gb,
                                        const ReciprocalLatticeVector<dim>& axis,
                                        const int& periodScaling):
        gb(gb),
        axis(axis),
        periodScaling(periodScaling)
    {
        auto basis= gb.bc.A.planeParallelLatticeBasis(gb.nA);
        ReciprocalLatticeDirection<3> planeOrthogonalToGb(  gb.bc.A.reciprocalLatticeDirection(shiftLatticeVector().cartesian() ));

        LatticeVector<3> axisCSL(gb.bc.csl.latticeDirection(axis.cartesian()).latticeVector());
        int numberOfCslPointsInGb= std::round(
                axisCSL.cartesian().cross(shiftLatticeVector().cartesian()).norm()/
                basis[1].cartesian().cross(basis[2].cartesian()).norm());

        IntScalarType alpha= basis[1].dot(planeOrthogonalToGb);
        IntScalarType beta= basis[2].dot(planeOrthogonalToGb);
        IntScalarType x,y;
        IntegerMath<IntScalarType>::extended_gcd(alpha, beta, x, y);
        LatticeVector<3> gbBasisAtom(x * basis[1].latticeVector() + y * basis[2].latticeVector());
        if(basis[0].dot(gb.nA)>0)
            basis[0]=-1*basis[0].latticeVector();

        const int stacking= gb.nA.stacking();
        for(int j=-stacking+1; j< stacking; ++j) {
            for (int i = 0; i < numberOfCslPointsInGb; ++i) {
                LatticeVector<3> atom((i+1)*gbBasisAtom);
                atom= atom + j*basis[0].latticeVector();
                Triplet index;
                index(0) =
                        IntegerMath<IntScalarType>::positive_modulo(atom.dot(planeOrthogonalToGb) , numberOfPlanesOrthogonalToGB());
                index(1) =
                        IntegerMath<IntScalarType>::positive_modulo(atom.dot(axis) , (ReciprocalLatticeDirection<3>(axis).stacking()));
                index(2) = j;
                refState.push_back(index);
            }
        }

    }

    template<int dim>
    int ReferenceState<dim>::numberOfPlanesOrthogonalToGB() const
    {
        int numberOfAPlanes= gb.bc.A.reciprocalLatticeDirection(shiftLatticeVector().cartesian()).stacking();
        int numberOfBPlanes= gb.bc.B.reciprocalLatticeDirection(shiftLatticeVector().cartesian()).stacking();
        return periodScaling*IntegerMath<int>::lcm(numberOfAPlanes,numberOfBPlanes);
    }

    template<int dim>
    LatticeVector<dim> ReferenceState<dim>::shiftLatticeVector() const {
        return periodScaling * gb.getPeriodVector(axis);
    }

    template class ReferenceState<3>;
} // gbLAB