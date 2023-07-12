#include <LatticeModule.h>
#include <TextFileParser.h>
#include <BiCrystal.h>

using namespace gbLAB;

int main()
{
    /*! [Types] */
    using VectorDimI = LatticeCore<3>::VectorDimI;
    using IntScalarType = LatticeCore<3>::IntScalarType;
    /*! [Types] */

    /*! [Lattice] */
    const auto A(TextFileParser("bicrystal_3d.txt").readMatrix<double,3,3>("A",true));
    Lattice<3> lattice(A);
    std::cout << "Lattice A = " << std::endl;
    std::cout << lattice.latticeBasis << std::endl;
    /*! [Lattice] */

    /*! [Axis] */
    const auto axis (TextFileParser("bicrystal_3d.txt").readMatrix<double,3,1>("axis",true));
    ReciprocalLatticeVector<3> rv(lattice.reciprocalLatticeDirection(axis).reciprocalLatticeVector());
    std::cout << "Cartesian coordinates of axis = " << std::endl;
    std::cout << rv.cartesian().transpose() << std::endl;
    /*! [Axis] */

    /*! [Generate bicrystal] */
    const auto& coincidentLattices= lattice.generateCoincidentLattices(rv,150,150);
    /*! [Generate bicrystal] */

    /*! [Misorientation] */
    for (const auto& rotation : coincidentLattices)
    {
        // Loop over misorientation angles
        std::cout << "###################################################" << std::endl;
        try
        {
            double theta= acos((rotation.trace()-1.0)/2.0)*180/M_PI;
            std::cout << "Misorientation angle = " << std::setprecision(20) << theta << "; ";
            BiCrystal<3> bc(lattice,Lattice<3>(lattice.latticeBasis,rotation),false);
            std::cout << "Sigma = " << std::setprecision(20) << bc.sigma << std::endl;
            std::cout << std::endl;
            std::cout << "Lattice B = " << std::endl;
            std::cout << std::setprecision(20) << rotation*lattice.latticeBasis << std::endl;
            std::cout << std::endl;
            std::cout << "Parallel CSL basis Cp= " << std::endl;
            std::cout << std::setprecision(20) << bc.csl.latticeBasis <<  std::endl;
            std::cout << std::endl;
            std::cout << "Parallel DSCL basis Dp = " << std::endl;
            std::cout << std::setprecision(20) << bc.dscl.latticeBasis <<  std::endl;
            std::cout << std::endl;

            // reduce DSCL basis vectors
            auto reducedDsclBasis= RLLL(bc.dscl.latticeBasis,0.75);
            auto U_Dscl= reducedDsclBasis.unimodularMatrix();

            std::cout << "Reduced DSCL basis vectors:" << std::endl;
            std::cout << "d1 = ";
            std::cout << std::setprecision(20) << reducedDsclBasis.reducedBasis().col(0).transpose() << std::endl;
            std::cout << "Integer coordinates of d1:";
            LatticeVector<3> d1(bc.dscl);
            d1 << U_Dscl.col(0).template cast<IntScalarType>();
            std::cout << std::setprecision(20) << d1.transpose() << std::endl;
            std::cout << std::endl;

            LatticeVector<3> d2(bc.dscl);
            std::cout << "d2 = ";
            std::cout << std::setprecision(20) << reducedDsclBasis.reducedBasis().col(1).transpose() << std::endl;
            std::cout << "Integer coordinates of d2:";
            d2 << U_Dscl.col(1).template cast<IntScalarType>();
            std::cout << std::setprecision(20) << d2.transpose() << std::endl;
            std::cout << std::endl;

            // shift vectors corresponding to d1 and d2
            LatticeVector<3> s1(bc.dscl), s2(bc.dscl);
            s1 << bc.LambdaA * d1;
            s2 << bc.LambdaA * d2;
            Lattice<3> reducedCsl(RLLL(bc.csl.latticeBasis,0.75).reducedBasis());
            std::cout << "Reduced shift vectors: " << std::endl;

            Eigen::Vector3d s1_coordinates_in_reduced_csl= reducedCsl.latticeBasis.inverse()*s1.cartesian();
            Eigen::Vector3d s2_coordinates_in_reduced_csl= reducedCsl.latticeBasis.inverse()*s2.cartesian();
            Eigen::Vector3d s1_coordinates_modulo= s1_coordinates_in_reduced_csl.array()-s1_coordinates_in_reduced_csl.array().round();
            Eigen::Vector3d s2_coordinates_modulo= s2_coordinates_in_reduced_csl.array()-s2_coordinates_in_reduced_csl.array().round();
            std::cout << "s1 = ";
            std::cout << std::setprecision(20) << (reducedCsl.latticeBasis * s1_coordinates_modulo).transpose() << std::endl;
            std::cout << "s2 = ";
            std::cout << std::setprecision(20) << (reducedCsl.latticeBasis * s2_coordinates_modulo).transpose() << std::endl;
            std::cout << std::endl;


            /*! [Misorientation] */
            /*! [Generate GBs] */
            auto gbSet(    bc.generateGrainBoundaries(bc.A.latticeDirection(rv.cartesian()),60) );
            /*! [Generate GBs] */
            /*! [Inclination] */
            int gbCount= 0;
            ReciprocalLatticeVector<3> refnA(bc.A);
            std::cout << "GBs of varying inclination (measured with respect to the first grain boundary" << std::endl;
            std::cout << "-----------------------------------------------------------------------------" << std::endl;
            for (const auto& gb : gbSet)
            {
                /*! [Inclination] */
                // Loop over inclinations
                /*! [Glide] */
                try
                {
                    LatticeVector<3> glide(rv.cross(gb.second.nA.reciprocalLatticeVector()).latticeVector());
                    LatticeVector<3> burgersVector(bc.getLatticeDirectionInD(glide).latticeVector());
                    LatticeVector<3> periodVector(bc.getLatticeDirectionInC(glide).latticeVector());

                    ReciprocalLatticeVector<3> rvInA= gb.second.nA.reciprocalLatticeVector();
                    ReciprocalLatticeVector<3> rvInCsl= bc.getReciprocalLatticeDirectionInC(rvInA).reciprocalLatticeVector();
                    /*! [Glide] */

                    /*! [Reference] */
                    if (gbCount==0) refnA= gb.second.nA.reciprocalLatticeVector();
                    /*! [Reference] */

                    /*! [Output] */
                    double epsilon=1e-8;
                    if (gbCount==0 || periodVector.cartesian().norm() < 100) {
                        double cosAngle= refnA.cartesian().normalized().dot(gb.second.nA.cartesian().normalized());
                        if (cosAngle-1>-epsilon) cosAngle= 1.0;
                        if (cosAngle+1<epsilon) cosAngle= -1.0;
                        std::cout << gbCount+1 << ") Inclination = " << std::setprecision(20) << acos(cosAngle)*180/M_PI << std::endl;
                        std::cout << "nA = " << gb.second.nA << std::endl;
                        std::cout << "nB = " << gb.second.nB << std::endl;
                        std::cout << "GB period = " << std::setprecision(20) << periodVector.cartesian().norm() << std::endl;
                        std::cout << "CSL plane distance (Height)= " << std::setprecision(20) << 1.0/rvInCsl.cartesian().norm() << std::endl;
                        std::cout << "Glide disconnection Burgers vector = " << std::setprecision(20) << burgersVector.cartesian().transpose() << "; norm = " << burgersVector.cartesian().norm() << std::endl;
                        std::cout << "Step height of glide disconnection = " << std::setprecision(20) << gb.second.stepHeightA(burgersVector) << std::endl;
                        std::cout << "Step height of non-glide disconnection 1= " << std::setprecision(20) << gb.second.stepHeightA(d1) << std::endl;
                        std::cout << "Step height of non-glide disconnection 2= " << std::setprecision(20) << gb.second.stepHeightA(d2) << std::endl;
                        std::cout << "-----------------------------------------------------------------------------" << std::endl;
                        gbCount++;
                    }
                    /*! [Output] */
                }
                catch(std::runtime_error& e)
                {
                    std::cout << e.what() << std::endl;
                    std::cout << "Moving onto the next inclination" << std::endl;
                }
            }
        }
        catch(std::runtime_error& e)
        {
            std::cout << e.what() << std::endl;
            std::cout << "Moving on the the next misorientation" << std::endl;
        }
    }
    return 0;
}
