#include <LatticeModule.h>
#include <TextFileParser.h>
#include <BiCrystal.h>
#include <MesoState.h>
#include <MesoStateEnsemble.h>

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
        /*! [SNF] */
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

        std::vector<std::pair<int,int>> dislocations;

        ReferenceState<3> rs(gb,rAxisA, 1);
        //MesoStateEnsemble<3> mesostates(rs,"mesostatesIn.txt");
        //MesoStateEnsemble<3> mesostates(rs,"builtMesostates.txt");
        MesoStateEnsemble<3> mesostatesForTraining(rs,100);

        std::cout << "Total number of training mesostates = " << mesostatesForTraining.size() << std::endl;
        MesoStateEnsemble<3> mesostates(mesostatesForTraining.build());
        mesostates.write("trainingStates.txt");
        //builtMesostates.write("builtMesostates.txt");
        //mesostates.write("mesostatesOut.txt");

        int mesostateIndex;
        std::ofstream energiesFile;
        energiesFile.open("energy.txt");
        energiesFile << "Mesostate" << std::setw(15) << "# CSL points" << std::setw(15) << "Elastic energy" << std::endl;
        mesostateIndex= 0;
        //for(const auto& ms : builtMesostates) {
        int numberOfInteractingPlanes= 15;
        for(const auto& ms : mesostates) {
            energiesFile << mesostateIndex << std::setw(15) << ms.getLocalStateCount(numberOfInteractingPlanes).transpose() << std::setw(15) << ms.elasticEnergy() << std::endl;
            mesostateIndex++;
        }
        energiesFile.close();
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
