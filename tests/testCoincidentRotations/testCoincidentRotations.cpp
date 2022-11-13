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
    /*! [Lattice] */

    /*! [Axis] */
    const auto axis (TextFileParser("bicrystal_3d.txt").readMatrix<IntScalarType,3,1>("axis",true));
    ReciprocalLatticeVector<3> rv(axis,lattice);
    ReciprocalLatticeDirection<3> rd(rv);
    std::cout << rv;
    /*! [Axis] */

    /*! [Test] */
    const auto& coincidentLattices= lattice.generateCoincidentLattices(rv);
    /*! [Test] */

    /*! [SNF] */
    for (const auto& rotatedLattice : coincidentLattices)
    {
        try
        {
            double theta= acos((rotatedLattice.second.F.trace()-1.0)/2.0)*180/M_PI;
            std::cout << "angle = " << theta << "; ";
            BiCrystal<3> bc(lattice,rotatedLattice.second,false);
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
