#include <TextFileParser.h>
#include <GbMesoStateEnsemble.h>
#include <MonteCarlo.h>
#include <LandauWangTP.h>

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


    // Sigma 29 [0-10](2 0 -5)
    VectorDimD axis(0,-1,0);
    double theta= 43.60282*M_PI/180;       // misorientation angle
    VectorDimD gbNormal(2,0,5);            // Miller indices
    int heightScaling= 1;
    int periodScaling= 1;
    int axisScaling= 1;
    double bScaling= 1.5;

    /*
    // Sigma 123 [110](-5 5 14)
    VectorDimD axis(1,1,0);
    double theta= 53.594515175286005615*M_PI/180;       // misorientation angle
    VectorDimD gbNormal(-5,5,14);                        // Miller indices
    int heightScaling= 2;
    int periodScaling= 1;
    int axisScaling= 1;
    double bScaling= 1.4;
     */

    /*
    //Sigma 3[1-10](112)
    VectorDimD axis(1,-1,0);
    double theta= 70.52878*M_PI/180;
    VectorDimD gbNormal(1,1,2);
    int heightScaling= 1;
    int periodScaling= 1;
    int axisScaling= 1;
    double bScaling= 2.0;
     */


    /*! [Lattice] */
    const auto A(TextFileParser("bicrystal_3d.txt").readMatrix<double,3,3>("A",true));
    Lattice<3> lattice(A);
    std::cout << "Lattice A = " << std::endl;
    std::cout << lattice.latticeBasis << std::endl;
    /*! [Lattice] */

    ReciprocalLatticeVector<3> rAxisGlobal(lattice.reciprocalLatticeDirection(axis).reciprocalLatticeVector());
    std::cout << "Cartesian coordinates of axis = " << std::endl;
    std::cout << rAxisGlobal.cartesian().transpose() << std::endl;

    try
    {
        // construct bicrystal
        Eigen::AngleAxis<double> halfRotation(theta/2,rAxisGlobal.cartesian().normalized());
        Lattice<3> latticeA(lattice.latticeBasis,halfRotation.matrix());
        Lattice<3> latticeB(lattice.latticeBasis,halfRotation.matrix().transpose());
        BiCrystal<3> bc(latticeA,latticeB,false);

        // construct GB
        gbNormal= halfRotation.matrix() * gbNormal;
        ReciprocalLatticeDirection<3> rd= latticeA.reciprocalLatticeDirection(gbNormal);
        Gb<3> gb(bc,rd);


        // Define the CSL box vectors
        ReciprocalLatticeVector<3> rAxisA(latticeA.reciprocalLatticeDirection(axis).reciprocalLatticeVector());
        std::cout << gb.getPeriodVector(rAxisA).cartesian().transpose() << std::endl;
        LatticeVector<3> axisA(gb.bc.A.latticeDirection(axis).latticeVector());
        LatticeVector<3> axisC(gb.bc.getLatticeDirectionInC(axisA).latticeVector());
        std::vector<LatticeVector<3>> cslVectors;
        cslVectors.push_back(heightScaling*gb.bc.csl.latticeDirection(gb.nA.cartesian()).latticeVector());
        cslVectors.push_back(periodScaling*gb.getPeriodVector(rAxisA));
        cslVectors.push_back(axisScaling*axisC);

        gb.box(cslVectors,1,1,"gb.txt");
        /*
         *  c11 = 1.0439923926128656 eV/angstrom^3
            c12 = 0.7750032094485771 eV/angstrom^3
            c44 = 0.4771664655306506 eV/angstrom^3/
            c11 = lambda + 2 mu;
            c12 = lambda; mu = (c11-c12)/2
            lambda/mu = 2 c12/(c11-c12)
            lattice constant = 3.615
         */
        // material parameter
        double a0= 3.615;
        double c11= 1.0439923926128656 * std::pow(a0,3);
        double c12= 0.7750032094485771 * std::pow(a0,3);
        GbMaterialTensors::lambda= c12;
        GbMaterialTensors::mu= (c11-c12)/2;

        GbMesoStateEnsemble<3> ensemble(gb, rAxisA, cslVectors, bScaling);

        /*
        auto mesostates= ensemble.collectMesoStates();
        for (const auto& mesostate : mesostates) {
            auto [density,energy]= mesostate.densityEnergy();
            std::cout << "density = " << density << "; energy = " << energy << std::endl;
        }
         */

        std::string potentialName= "Cu_mishin1.eam.alloy";
        std::string lmpLocation= "/Users/Nikhil/Documents/Academic/Software/lammps-15May15/src/lmp_serial";
        //LandauWangTP<XTuplet,GbMesoState<3>> landauWang(1.51, 20,30,0.78,0.96,9);
        LandauWangTP<XTuplet,GbMesoState<3>> landauWang({1.51, 20,30}, lmpLocation, potentialName);

        MonteCarlo<XTuplet, GbMesoState<3>, GbMesoStateEnsemble<3>, LandauWangTP<XTuplet,GbMesoState<3>>> mc(ensemble, landauWang);
        for (int i=0; i<50; ++i) {
            //const auto& constraintsMesostateMap= mc.evolve(1000,20000,"ms");
            mc.evolve(10000);
            std::ofstream outputFileHandle;
            outputFileHandle.open("theta"+std::to_string(i)+".txt");
            outputFileHandle << landauWang.theta;
        }


    }
    catch(std::runtime_error& e)
    {
        //std::cout << e.what() << std::endl;
    }
    return 0;
}
