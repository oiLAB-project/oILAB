#include <LatticeModule.h>
#include <TextFileParser.h>
#include <BiCrystal.h>

using namespace gbLAB;

int main()
{
    using VectorDimI = LatticeCore<3>::VectorDimI;
    using IntScalarType = LatticeCore<3>::IntScalarType;
    const auto A(TextFileParser("bicrystal_3d.txt").readMatrix<double,3,3>("A",true));
    const auto axis (TextFileParser("bicrystal_3d.txt").readMatrix<IntScalarType,3,1>("axis",true));
    Lattice<3> lattice(A);

    ReciprocalLatticeVector<3> rv(axis,lattice);
    ReciprocalLatticeDirection<3> rd(rv);
    const auto& coincidentLattices= lattice.generateCoincidentLattices(rv);
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
    return 0;
}
