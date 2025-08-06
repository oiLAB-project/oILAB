#include <LatticeModule.h>
#include <TextFileParser.h>

using namespace gbLAB;

int main()
{
    /*! [Types] */
    using VectorDimI = LatticeCore<3>::VectorDimI;
    using VectorDimD = LatticeCore<3>::VectorDimD;
    using IntScalarType = LatticeCore<3>::IntScalarType;
    /*! [Types] */

    /*! [Lattice] */
    const auto A(TextFileParser("bicrystal_3d.txt").readMatrix<double,3,3>("A",true));
    Lattice<3> lattice(A);
    /*! [Lattice] */

    /*! [Axis] */
    const auto axis (TextFileParser("bicrystal_3d.txt").readMatrix<double,3,1>("axis",true));
    ReciprocalLatticeVector<3> rv(axis,lattice);
    std::cout << "Miller indices of the tilt axis plane = " << std::endl;
    std::cout << rv << std::endl;
    /*! [Axis] */

    /*! [Test] */
    const auto& coincidentRotations= lattice.generateCoincidentLattices(rv);
    /*! [Test] */

    /*! [SNF] */
    for (const auto& rotation : coincidentRotations)
    {
        try
        {
            double theta= acos((rotation.trace()-1.0)/2.0)*180/M_PI;
            std::cout << "angle = " << theta << "; ";
            BiCrystal<3> bc(lattice,Lattice<3>(lattice.latticeBasis,rotation),false);
            std::cout << "Sigma = " << bc.sigma << std::endl;
        }
        catch(std::runtime_error& e)
        {
            std::cout << e.what() << "Moving on ..." << std::endl;
        }
    }
    /*! [SNF] */
    return 0;
}
