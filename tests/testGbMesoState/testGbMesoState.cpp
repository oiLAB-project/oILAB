#include <TextFileParser.h>
#include <GbMesoStateEnsemble.h>

using namespace gbLAB;

int main()
{
    /*! [Types] */
    using VectorDimI = LatticeCore<3>::VectorDimI;
    using VectorDimD = LatticeCore<3>::VectorDimD;
    using Vector2d= Eigen::Vector2d;
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
    ReciprocalLatticeVector<3> rAxisGlobal(lattice.reciprocalLatticeDirection(axis).reciprocalLatticeVector());
    std::cout << "Cartesian coordinates of axis = " << std::endl;
    std::cout << rAxisGlobal.cartesian().transpose() << std::endl;
    /*! [Axis] */

    try
    {
        double theta= 53.594515175286005615*M_PI/180;
        //double theta= 50.478803641357835374*M_PI/180;
        Eigen::AngleAxis<double> halfRotation(theta/2,rAxisGlobal.cartesian().normalized());
        Lattice<3> latticeA(lattice.latticeBasis,halfRotation.matrix());
        Lattice<3> latticeB(lattice.latticeBasis,halfRotation.matrix().transpose());

        BiCrystal<3> bc(latticeA,latticeB,false);

        // Specify GB normal
        VectorDimD gbNormal;
        gbNormal << 0.0, 0.0, 1.0;
        ReciprocalLatticeDirection<3> rd= latticeA.reciprocalLatticeDirection(gbNormal);

        Gb<3> gb(bc,rd);
        ReciprocalLatticeVector<3> rAxisA(latticeA.reciprocalLatticeDirection(axis).reciprocalLatticeVector());
        /*
         *  c11 = 1.0439923926128656 eV/angstrom^3
            c12 = 0.7750032094485771 eV/angstrom^3
            c44 = 0.4771664655306506 eV/angstrom^3/
            c11 = lambda + 2 mu;
            c12 = lambda; mu = (c11-c12)/2
            lambda/mu = 2 c12/(c11-c12)
            lattice constant = 3.615
         */

        //GbMaterialTensors::lambda= 2 * 0.7750032094485771 / (1.0439923926128656 - 0.7750032094485771);
        //GbContinuum<3>::lambda= 0.5;
        //GbMaterialTensors::mu= 1;

        double a0= 3.615;
        double c11= 1.0439923926128656 * std::pow(a0,3);
        double c12= 0.7750032094485771 * std::pow(a0,3);
        GbMaterialTensors::lambda= c12;
        GbMaterialTensors::mu= (c11-c12)/2;
        // The last argument 0.3 is optional, and it defaults to 0.3. Increasing this number will increase the number of mesostates
        GbMesoStateEnsemble<3> ensemble(gb, rAxisA, 0.3,{1,1,2});
        ensemble.collectMesoStates("ms");

    }
    catch(std::runtime_error& e)
    {
        std::cout << e.what() << std::endl;
        std::cout << "Moving on the the next misorientation" << std::endl;
        std::cout << "-----------------------------------------------------------------------------"
                  << std::endl;
    }
    return 0;
}
