//
// Created by Nikhil Chandra Admal on 12/11/23.
//

#ifndef OILAB_MESOSTATE_H
#define OILAB_MESOSTATE_H

#include <Gb.h>

namespace gbLAB {
    template<int dim>
    class MesoState {
        using IntScalarType= typename LatticeCore<dim>::IntScalarType;
        using VectorDimI= typename LatticeCore<dim>::VectorDimI;
        using Matrix= typename Eigen::Matrix<IntScalarType,Eigen::Dynamic,Eigen::Dynamic>;

    private:
        static LatticeVector<dim> getPeriodVector();
        static void getRefState();

    public:
        const Gb<dim>& gb;
        const ReciprocalLatticeDirection<dim>& axis;
        LatticeVector<dim> periodVector;
        Matrix refState, currentState;

        explicit MesoState(const Gb<dim>& gb, const LatticeDirection<dim>& axis) :
            gb(gb),
            axis(axis),
            periodVector(getPeriodVector()),
            refState(getRefState()),
            currentState(refState)
        {
                LatticeVector<3> glideA(axis.cross(gb.nA.reciprocalLatticeVector()).latticeVector());
                LatticeVector<3> burgersVectorA(gb.bc.getLatticeDirectionInD(glideA).latticeVector());
                LatticeVector<3> periodVectorA(gb.bc.getLatticeDirectionInC(glideA).latticeVector());

                ReciprocalLatticeVector<3> gbNormalA = gb.nA.reciprocalLatticeVector();
                ReciprocalLatticeVector<3> gbNormalCsl= gb.bc.getReciprocalLatticeDirectionInC(gbNormalA).reciprocalLatticeVector();

                /*! [Output] */
                double epsilon = 1e-8;
                if (periodVectorA.cartesian().norm() < 50) {
                    std::cout << "nA = " << gb.nA << std::endl;
                    std::cout << "nB = " << gb.nB << std::endl;
                    std::cout << "GB period = " << std::setprecision(20) << periodVectorA.cartesian().norm()
                              << std::endl;
                    std::cout << "CSL plane distance (Height)= " << std::setprecision(20)
                              << 1.0 / gbNormalCsl.cartesian().norm() << std::endl;

                    // planes perpendicular to the GB
                    ReciprocalLatticeDirection<3> planeOrthogonalToGb(  gb.bc.A.reciprocalLatticeDirection(burgersVectorA.cartesian() ));
                    int numberOfPlanesOrthogonalToGB= planeOrthogonalToGb.stacking();
                    std::cout << "Number of lattice planes perpendicular to the GB = "
                              << numberOfPlanesOrthogonalToGB << std::endl;
                    std::cout << "Number of lattice planes parallel to the GB = "
                              << gb.nA.stacking() << std::endl;



                    LatticeVector<3> axisCSL(gb.bc.csl.latticeDirection(axis).latticeVector());
                    auto basis= gb.bc.A.planeParallelLatticeBasis(gb.nA);
                    int numberOfCslPoints= std::round(
                            axisCSL.cartesian().cross(periodVectorA.cartesian()).norm()/
                            basis[1].cartesian().cross(basis[2].cartesian()).norm());
                    std::cout << "Number of CSL points = " << numberOfCslPoints << std::endl;

                    IntScalarType alpha= basis[1].dot(planeOrthogonalToGb);
                    IntScalarType beta= basis[2].dot(planeOrthogonalToGb);
                    IntScalarType x,y;
                    IntegerMath<IntScalarType>::extended_gcd(alpha, beta, x, y);
                    LatticeVector<3> gbBasisAtom(x * basis[1].latticeVector() + y * basis[2].latticeVector());
                    std::vector<VectorDimI> basisAtomIndex;

                    for(int j=0; j< gb.nA.stacking(); ++j) {
                        for (int i = 0; i < numberOfCslPoints; ++i) {
                            LatticeVector<3> atom((i+1)*gbBasisAtom);
                            atom= atom + j*basis[0].latticeVector();
                            VectorDimI index;
                            index(0) =
                                    IntegerMath<IntScalarType>::positive_modulo(atom.dot(planeOrthogonalToGb) , numberOfPlanesOrthogonalToGB);
                            index(1) =
                                    IntegerMath<IntScalarType>::positive_modulo(atom.dot(axis) , (ReciprocalLatticeDirection<3>(axis).stacking()));
                            index(2) = j;
                            basisAtomIndex.push_back(index);
                        }
                    }

                    std::map<int,int> mesostate;
                    for(auto index : basisAtomIndex)
                    {
                        mesostate.insert(std::make_pair<int,int>(index(0),index(2)));
                    }
                    for(auto pair : mesostate)
                        std::cout << pair.first << "  " << pair.second << std::endl;

                    // mesostate generator
                    std::vector<std::map<int,int>> mesostates;
                    mesostates.push_back(mesostate);


                    std::cout << "-----------------------------------------------------------------------------"
                              << std::endl;
            }
            /*! [Output] */



        }

        void transition();
        double surfaceEnergy();
        double elasticEnergy();
    };
}
#endif //OILAB_MESOSTATE_H
