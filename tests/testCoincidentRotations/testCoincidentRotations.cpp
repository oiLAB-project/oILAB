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
    auto coincidentLattices= lattice.generateCoincidentLattices(rv);
    for (const auto& rotatedLattice : coincidentLattices)
    {
        BiCrystal<3> bc(lattice,rotatedLattice,false);
        std::cout << bc.sigma << std::endl;
    }
    return 0;
}
