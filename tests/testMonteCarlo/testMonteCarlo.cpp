#include <LatticeModule.h>
#include <TextFileParser.h>
#include <BiCrystal.h>
#include <MesoState.h>
#include <MesoStateEnsemble.h>
#include <randomInteger.h>

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
        //rs.planeEnergies = (Eigen::VectorXd(15) << 3.05292806, 1027.29054285,  340.04377908,  181.13969142,
        //                    99.22400022,   27.9608597 ,  -31.67405067,  -17.28199054,
        //                    -25.65920418,  -12.84802797,  -12.28947694,   13.03675517,
        //                    -34.58309914,  -76.66654383,  -34.78187716).finished();
        rs.planeEnergies = (Eigen::VectorXd(10) << -10, 0,  0,  0,
                0, 0,  0, 0, 0, 0).finished();
        //rs.planeEnergies = (Eigen::VectorXd(10) << -2.26570104, 1020.90169538,  338.74619593,  120.68298265,
        //        53.68030732,    2.96439155,  -46.69213003,  -24.83383622,
        //        -24.35190955,   -8.79339053).finished();
        MesoStateEnsemble<3> mesostates(rs);

        MesoState<3> ms(rs, 1, 5);
        //mesostates.insert(ms);
        mesostates.push_back(ms);
        const int maxIterations= 1000000;
        const int maxNumberOfDipoles= 12;


        int mesostateIndex;
        std::ofstream energiesFile;
        energiesFile.open("energy.txt");
        energiesFile << "Mesostate" << std::setw(15) << "# CSL points" << std::setw(15) << "Elastic energy" << std::endl;
        int numberOfInteractingPlanes= rs.planeEnergies.size();
        mesostateIndex= 0;


        std::random_device r_d;
        std::mt19937 gen(r_d());
        for(int i=0; i< maxIterations; ++i) {
            MesoState<3> temp(ms);
            Triplet t;
            bool insert;
            if (ms.defectsIndices.size() == 0) {
                t = temp.insertRandomDislocation();
                insert=true;
            }
            else if (ms.defectsIndices.size() == maxNumberOfDipoles) {
                t = temp.removeRandomDislocation();
                insert = false;
            }
            else {
                if (random<int>(-100, 100) <= 0) {
                    t = temp.removeRandomDislocation();
                    insert = false;
                }
                else {
                    t = temp.insertRandomDislocation();
                    insert = true;
                }
            }

            double delta = temp.energy()-ms.energy();
            double temperature= 2.5;
            double probability = std::min(1.0,exp(-delta/temperature));

            if (random<double>(0.0,1.0)<=probability) {
            //if (delta<=0){
                if (insert)
                    ms.insertDislocation(t);
                else
                    ms.removeDislocation(t);
                mesostates.push_back(ms);
                energiesFile << mesostateIndex << std::setw(15) <<
                             ms.getLocalStateCount(numberOfInteractingPlanes).transpose() << std::setw(15) <<
                             ms.energy() << std::setw(15) << std::endl;
                mesostateIndex++;
            }
        }

        energiesFile.close();
        std::cout << "Total number of mesostates = " << mesostates.size() << std::endl;
        mesostates.write("mesostates.txt");

        /*
        //for(const auto& ms : builtMesostates) {
        int numberOfInteractingPlanes= rs.planeEnergies.size();
        for(const auto& ms : mesostates) {
            energiesFile << mesostateIndex << std::setw(15) <<
            ms.getLocalStateCount(numberOfInteractingPlanes).transpose() << std::setw(15) <<
            ms.energy() << std::setw(15) << std::endl;
            //ms.elasticEnergy() << std::endl;
            mesostateIndex++;
        }
        energiesFile.close();
         */


        /*
        MesoStateEnsemble<3> builtMesoStates(rs);
        double temp;
        temp= builtMesoStates[0].energy();
        for (const auto& ms : mesostates) {
            //if (temp-ms.energy() > 0.5*abs(temp)) {
                temp = ms.energy();
                builtMesoStates.push_back(ms);
            //}
        }
        //auto builtMesoStates(mesostates.build());
        std::cout << "Size of built mesostates  = " << builtMesoStates.size() << std::endl;
        builtMesoStates.build();
        builtMesoStates.write("builtMesostates.txt");
         */
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
