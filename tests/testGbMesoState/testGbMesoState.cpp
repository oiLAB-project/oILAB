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
    int heightScaling= 2;
    int periodScaling= 2;
    int axisScaling= 1;
    double bScaling= 2.0;
     */


    /*! [Lattice] */
    double a0= 3.615000084042549;
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

	std::deque<XTuplet> constraintsEnsemble(ensemble.enumerateConstraints((const GbShifts<3>&) ensemble));
	std::cout << "Size of the ensemble = " << constraintsEnsemble.size() << std::endl;

        std::ofstream out_file;

        #pragma omp parallel for num_threads(20) private(out_file)
        for (size_t i = 0; i < constraintsEnsemble.size(); ++i) {
            int thread_id = omp_get_thread_num();

            // Create a filename for each thread
            std::string filename = "output_thread_" + std::to_string(thread_id) + ".txt";

            // Open the file once per thread (if not already open)
            if (!out_file.is_open()) 
                out_file.open(filename);

            const auto& it = std::next(constraintsEnsemble.begin(), i);
            try {
                const auto& mesostate= ensemble.constructMesoState(*it);
                const auto data= mesostate.densityEnergy();
	        if (out_file.is_open())
                    out_file << *it << "  " << it->density() << "  " << data.second << std::endl;
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
        //std::cout << e.what() << std::endl;
    }
    return 0;
}
