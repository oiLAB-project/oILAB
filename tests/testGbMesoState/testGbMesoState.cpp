#include <TextFileParser.h>
#include <GbMesoStateEnsemble.h>

using namespace gbLAB;

// Hunter GB Details
//
// Sigma 3 [110]; theta = 70.52878
// GB normal: {1 - 1 2}

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
        double theta= 70.52878*M_PI/180;
        //double theta= 53.594515175286005615*M_PI/180;
        //double theta= 43.60282*M_PI/180;
        Eigen::AngleAxis<double> halfRotation(theta/2,rAxisGlobal.cartesian().normalized());
        Lattice<3> latticeA(lattice.latticeBasis,halfRotation.matrix());
        Lattice<3> latticeB(lattice.latticeBasis,halfRotation.matrix().transpose());

        BiCrystal<3> bc(latticeA,latticeB,false);

        // Specify GB normal
        VectorDimD gbNormal;
        //gbNormal << 0.0, 0.0, 1.0;
        gbNormal << 1.0, 1.0, 2.0;
        gbNormal= halfRotation.matrix() * gbNormal;
        ReciprocalLatticeDirection<3> rd= latticeA.reciprocalLatticeDirection(gbNormal);

        Gb<3> gb(bc,rd);
        ReciprocalLatticeVector<3> rAxisA(latticeA.reciprocalLatticeDirection(axis).reciprocalLatticeVector());
        LatticeVector<3> axisA(gb.bc.A.latticeDirection(axis).latticeVector());
        LatticeVector<3> axisC(gb.bc.getLatticeDirectionInC(axisA).latticeVector());

        std::vector<LatticeVector<3>> cslVectors;
        cslVectors.push_back(15*gb.bc.csl.latticeDirection(gb.nA.cartesian()).latticeVector());
        cslVectors.push_back(2*gb.getPeriodVector(rAxisA));
        cslVectors.push_back(1*axisC);

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
        GbMesoStateEnsemble<3> ensemble(gb, rAxisA, cslVectors, 2.2);
        ensemble.collectMesoStates("ms");

    }
    catch(std::runtime_error& e)
    {
        std::cout << e.what() << std::endl;
    }
    return 0;
}
