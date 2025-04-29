#include <TextFileParser.h>
#include <GbMesoStateEnsemble.h>
#include <omp.h>

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

    /*
    // Sigma 37[111], misorientation = 50.569992092103568382
    VectorDimD axis(1,1,1);
    double theta = 50.569992092103568382*M_PI/180;
    VectorDimD gbNormal(10,1,-11);
    int heightScaling= 1;
    int periodScaling= 1;
    int axisScaling= 1;
    double bScaling= 0.75;
     */


    /*
    // Sigma 29 [0-10](2 0 -5)
    VectorDimD axis(0,-1,0);
    double theta= 43.60282*M_PI/180;       // misorientation angle
    VectorDimD gbNormal(2,0,5);            // Miller indices
    int heightScaling= 1;
    int periodScaling= 1;
    int axisScaling= 1;
    double bScaling= 0.75;
     */

    // Sigma 123 [110](-5 5 14)
    VectorDimD axis(1,1,0);
    double theta= 53.594515175286005615*M_PI/180;       // misorientation angle
    VectorDimD gbNormal(-5,5,14);                        // Miller indices
    int heightScaling= 2;
    int periodScaling= 1;
    int axisScaling= 1;
    double bScaling= 0.6;

    //Sigma 3
    /*
    VectorDimD axis(1,-1,0);
    double theta= 70.52878*M_PI/180; //Sigma 3[1-10](112) incoherent
    VectorDimD gbNormal(1,1,2);
    //double theta= 109.47122*M_PI/180; //Sigma 3[1-10](111) coherent
    //VectorDimD gbNormal(1,1,1);
    int heightScaling= 4;
    int periodScaling= 1;
    int axisScaling= 1;
    double bScaling= 0.99;
     */

    /*
    //Sigma 11
    VectorDimD axis(1,-1,0);
    double theta= 50.478803641357835374*M_PI/180; //Sigma 3[1-10](112) incoherent
    VectorDimD gbNormal(1,1,3);
    int heightScaling= 2;
    int periodScaling= 1;
    int axisScaling= 1;
    double bScaling= 0.99;
    */


    /*! [Lattice] */
    double strain=0.0; // 0.01 or 0.02
    double a0= (1.0+strain)*3.615000084042549;
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

    try
    {
        // construct bicrystal
        Eigen::AngleAxis<double> halfRotation(theta/2,rAxisGlobal.cartesian().normalized());
        Lattice<3> latticeA(lattice.latticeBasis,halfRotation.matrix());
        Lattice<3> latticeB(lattice.latticeBasis,halfRotation.matrix().transpose());
        BiCrystal<3> bc(latticeA,latticeB,false);
        std::cout << "Sigma = " << bc.sigma << std::endl;

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
        gb.box(cslVectors,1,1,"gb.txt",true);
        bc.box(cslVectors,1,1,"bcOriented.txt",true);


        // material parameter
        // source: https://openkim.org/id/EAM_Dynamo_MishinMehlPapaconstantopoulos_2001_Cu__MO_346334655118_005
        double c11= 169.9281940954852/160.2176621;
        double c12= 122.65063014404001/160.2176621;
        GbMaterialTensors::lambda= c12;
        GbMaterialTensors::mu= (c11-c12)/2;

        GbMesoStateEnsemble<3> ensemble(gb, rAxisA, cslVectors, bScaling);
        std::deque<XTuplet> constraintsEnsemble(ensemble.enumerateConstraints((const GbShifts<3>&) ensemble));
        std::cout << "Size of the ensemble = " << constraintsEnsemble.size() << std::endl;

        std::ofstream out_file;

        #pragma omp parallel for num_threads(1) private(out_file)
        for (size_t i = 0; i < constraintsEnsemble.size(); ++i) {
            int thread_id = omp_get_thread_num();

            // Create a filename for each thread
            std::string filename= "output_thread_" + std::to_string(thread_id) + ".txt";

            // Open the file once per thread (if not already open)
            if (!out_file.is_open()) 
                out_file.open(filename);

            const auto& it = std::next(constraintsEnsemble.begin(), i);
            try {
                const auto& mesostate= ensemble.constructMesoState(*it);
                //mesostate.box(std::to_string(i)+ ".txt");
                std::string potentialName= "Cu_mishin1.eam.alloy";
                std::string lmpLocation= "/Users/Nikhil/Documents/Academic/Software/lammps-15May15/src/lmp_serial";

                const auto data= mesostate.densityEnergy(lmpLocation, potentialName, true, {10,10,10});
                PeriodicFunction<double,3> rho= std::get<2>(data) ;
                // the density values can be read using rho.values

                if (out_file.is_open())
                        //out_file << *it << "  " << it->density() << "  " << data.second << std::endl;
                        out_file << *it << "  " << std::get<0>(data)<< "  " << std::get<1>(data) << std::endl;
                else
                    std::cerr << "Failed to open file " << filename << std::endl;
            }
            catch(std::runtime_error& e)
            {
                std::cout << e.what() << std::endl;
            }
        }
    }
    catch(std::runtime_error& e)
    {
        std::cout << e.what() << std::endl;
    }
    return 0;
}
