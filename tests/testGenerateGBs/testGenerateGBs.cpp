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
    const auto& coincidentLattices= lattice.generateCoincidentLattices(rv);
    /*! [Generate bicrystal] */

    /*! [Misorientation] */
    for (const auto& rotatedLattice : coincidentLattices)
    {
        // Loop over misorientation angles
        std::cout << "###################################################" << std::endl;
        try
        {
            double theta= acos((rotatedLattice.second.F.trace()-1.0)/2.0)*180/M_PI;
            std::cout << "Misorientation angle = " << theta << "; ";
            BiCrystal<3> bc(lattice,rotatedLattice.second,false);
            std::cout << "Sigma = " << bc.sigma << std::endl;
            std::cout << "Lattice B = " << std::endl;
            std::cout << rotatedLattice.second.latticeBasis << std::endl;
            /*! [Misorientation] */
            /*! [Generate GBs] */
            auto gbSet(bc.generateGrainBoundaries(bc.A.latticeDirection(rv.cartesian())));
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
                        std::cout << gbCount+1 << ") Inclination = " << acos(cosAngle)*180/M_PI << std::endl;
                        std::cout << "nA = " << gb.second.nA << std::endl;
                        std::cout << "nB = " << gb.second.nB << std::endl;
                        std::cout << "GB period = " << periodVector.cartesian().norm() << std::endl;
                        std::cout << "Burgers vector = " << burgersVector.cartesian().transpose() << "; norm = " << burgersVector.cartesian().norm() << std::endl;
                        std::cout << "Step height of glide dislocation = " << gb.second.stepHeightA(burgersVector)
                                  << std::endl;
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
