#include <TextFileParser.h>
#include <GbMesoStateEnsemble.h>
#include <MonteCarlo.h>
#include <CanonicalTP.h>
#include <omp.h>

using namespace gbLAB;

void runMonteCarlo(const double& a0,
                   const double& temperature,
                   const int& iterations,
                   XTuplet startingState,
                   const std::string filename="")
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
    double bScaling= 2.0;

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
    Eigen::Matrix3d A;
    A << 0.0, 0.5, 0.5,
            0.5, 0.0, 0.5,
            0.5, 0.5, 0.0;
    A= a0*A;
    Lattice<3> lattice(A);
    std::cout << "Lattice A = " << std::endl;
    std::cout << lattice.latticeBasis << std::endl;
    /*! [Lattice] */

    ReciprocalLatticeVector<3> rAxisGlobal(lattice.reciprocalLatticeDirection(axis).reciprocalLatticeVector());
    std::cout << "Cartesian coordinates of axis = " << std::endl;
    std::cout << rAxisGlobal.cartesian().transpose() << std::endl;

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

    // material parameter
    // source: https://openkim.org/id/EAM_Dynamo_MishinMehlPapaconstantopoulos_2001_Cu__MO_346334655118_005
    double c11= 169.9281940954852/160.2176621;
    double c12= 122.65063014404001/160.2176621;
    GbMaterialTensors::lambda= c12;
    GbMaterialTensors::mu= (c11-c12)/2;

    GbMesoStateEnsemble<3> ensemble(gb, rAxisA, cslVectors, bScaling);
    const double kb = 8.61733e-5; // in eV/K
    CanonicalTP<XTuplet, GbMesoState<3>> canonicalTP(kb*temperature,filename);
    if (startingState.size()!=0) {
        MonteCarlo<XTuplet, GbMesoState < 3>, GbMesoStateEnsemble < 3 >, CanonicalTP<XTuplet, GbMesoState < 3>>> mc(
                ensemble, canonicalTP, startingState);
        std::cout << "Starting MC with state = " << startingState << std::endl;
        mc.evolve(iterations);
    }
    else {
        MonteCarlo<XTuplet, GbMesoState < 3>, GbMesoStateEnsemble < 3 >, CanonicalTP<XTuplet, GbMesoState < 3>>> mc(
                ensemble, canonicalTP,ensemble.initializeState());
        mc.evolve(iterations);
    }
    GbContinuum<3>::reset();
}

std::pair<double,XTuplet> getLowestEnergyState(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
    }
    std::string line;
    std::vector<std::string> row_with_min_energy;
    double min_value = std::numeric_limits<double>::infinity(); // Initialize with a very high value

// Read file line by line
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::vector<std::string> row;
        std::string value;

// Split the line into individual columns
        while (ss >> value) {
            row.push_back(value);
        }

// Check if the row is not empty and has at least one value
        if (!row.empty()) {
            try {
// Get the last column value (last element in the row)
                double current_value = std::stod(row.back());

// Check if this value is the lowest we have encountered
                if (current_value < min_value) {
                    min_value = current_value;
                    row_with_min_energy = row; // Store the entire row with the minimum value
                }
            } catch (const std::invalid_argument &e) {
// Handle the case where conversion to double fails
                std::cerr << "Invalid number in the file: " << e.what() << std::endl;
            }
        }
    }
    file.close();

    double lowestEnergy;
    std::vector<int> lowestEnergyStateVector;
    if (!row_with_min_energy.empty()) {
        lowestEnergy= std::stod(row_with_min_energy.back());
        for (auto it = row_with_min_energy.begin(); it != std::prev(std::prev(row_with_min_energy.end())); ++it) {
            lowestEnergyStateVector.push_back(std::stoi(*it));
        }
    } else {
        std::cerr << "No valid data found in the file." << std::endl;
    }

    XTuplet lowestEnergyState(lowestEnergyStateVector.size());
    for(int i=0; i<lowestEnergyStateVector.size(); ++i)
        lowestEnergyState(i)= lowestEnergyStateVector[i];

    return std::make_pair(lowestEnergy,lowestEnergyState);
}

int main()
{

    // run a high temperature monte carlo
    double a0= 3.615000084042549;
    double temperature= 3000;
    int iterations= 1;
    std::string filename = "lowestEnergy";
    runMonteCarlo(a0,temperature,iterations,XTuplet(0),filename);


    // get lowest energy state
    auto pair= getLowestEnergyState(filename);
    std::cout << "Lowest energy = " << pair.first << std::endl;
    std::cout << "Lowest energy state = " << pair.second << std::endl;
    XTuplet lowestEnergyState(pair.second);
    lowestEnergyState << 1, 1, 1, 1, 1, 1, 1, 2, 0, 0, 2, 0, 2, 2;
    //lowestEnergyState << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    std::cout << lowestEnergyState << std::endl;


    double a0ref= a0;
    iterations = 1;
#pragma omp parallel for num_threads(3) collapse(2)
    for(int i=1; i<=3; ++i)
    {
        for(int j=0; j<=10; ++j)
        {
            temperature = 300 * i;
	        double strain= (double)j/500.0;
            a0= (1.0+strain)*a0ref;
            runMonteCarlo(a0, temperature, iterations, lowestEnergyState, "observablesT" + std::to_string(temperature) + "_a" + std::to_string(strain));
        }
    }
    std::cout << "Done" << std::endl;
    return 0;
}
