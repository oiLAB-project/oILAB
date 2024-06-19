//
// Created by Nikhil Chandra Admal on 2/4/24.
//
#include <GbShifts.h>
#include <randomInteger.h>

namespace gbLAB
{
    template<int dim>
    GbShifts<dim>::GbShifts(const Gb<dim>& gb,
                            const ReciprocalLatticeVector<dim>& axis,
                            const double& bhalfMax):
            gb(gb),
            axis(axis),
            gbCslVectors(getGbCslVectors(gb,axis))
    {
        // assert that axis should belong to one of the four reciprocal lattices
        auto nC= gb.bc.getReciprocalLatticeDirectionInC(gb.nB.reciprocalLatticeVector());
        ReciprocalLatticeVector<dim> raxisD(gb.bc.getReciprocalLatticeDirectionInD(axis).reciprocalLatticeVector());
        ReciprocalLatticeDirection<dim> raxisT(gb.getReciprocalLatticeDirectionInT(raxisD));

        // CSL sublattice formed by 1) period vector, 2) axis, and 3) out-plane csl vector
        auto gbPlaneParallelCslBasis= gb.bc.csl.planeParallelLatticeBasis(nC,true);
        LatticeVector<dim> axisA(gb.bc.A.latticeDirection(axis.cartesian()).latticeVector());
        LatticeVector<dim> axisC(gb.bc.getLatticeDirectionInC(axisA).latticeVector());

        std::vector<LatticeVector<dim>> cslSubLatticeVectors;
        cslSubLatticeVectors.push_back(gbPlaneParallelCslBasis[0].latticeVector());
        cslSubLatticeVectors.push_back(gbCslVectors[0]);
        cslSubLatticeVectors.push_back(gbCslVectors[1]);
        //cslSubLatticeVectors.push_back(gb.getPeriodVector(axis));
        //cslSubLatticeVectors.push_back(axisC);

        // T Sublattice
        std::vector<LatticeVector<dim>> latticeVectorsT;
        int numPlanes= floor(bhalfMax/gb.nA.planeSpacing());
        std::cout << "GBShifts: Number of planes explored =  " << numPlanes << std::endl;
        latticeVectorsT.push_back(numPlanes * gb.T.latticeDirection(gb.nA.cartesian()).latticeVector());
        latticeVectorsT.push_back(gb.T.latticeDirection(gb.getPeriodVector(axis).cartesian()).latticeVector());
        latticeVectorsT.push_back(gb.T.latticeDirection(axisC.cartesian()).latticeVector());

        auto points= gb.T.box(latticeVectorsT);

        VectorDimD shiftT, shiftC;
        shiftT << -0.5, -0.5, -0.5;
        shiftC << -0.5,  0.0, -0.5;
        for(auto& point : points) {
            LatticeVector<dim>::modulo(point, latticeVectorsT, shiftT);
            auto cslShift = LatticeVector<dim>((gb.bc.LambdaA * gb.basisT * point).eval(), gb.bc.dscl);
            VectorDimD cslShiftCentered = cslShift.cartesian() - point.cartesian() / 2;
            LatticeVector<dim>::modulo(cslShiftCentered, cslSubLatticeVectors, shiftC);
            bShiftPairs.push_back(std::make_pair(point, cslShiftCentered));
        }
    }

    template<int dim>
    std::vector<LatticeVector<dim>> GbShifts<dim>::getGbCslVectors(const Gb<dim>& gb, const ReciprocalLatticeVector<dim>& axis)
    {
        std::vector<LatticeVector<dim>> output;
        LatticeVector<dim> axisA(gb.bc.A.latticeDirection(axis.cartesian()).latticeVector());
        LatticeVector<dim> axisC(gb.bc.getLatticeDirectionInC(axisA).latticeVector());
        for(int i=0; i<dim-1; ++i) {
            if (i==0) output.push_back(gb.getPeriodVector(axis));
            if (i==1) output.push_back(axisC);
        }
        return output;
    }

    //template class GbShifts<2>;
    template class GbShifts<3>;

}
