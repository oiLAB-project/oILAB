#include <LatticeModule.h>
#include <TextFileParser.h>
#include <BiCrystal.h>
#include <chrono>

using namespace gbLAB;

int main()
{
	auto start = std::chrono::system_clock::now();
	std::time_t startTime = std::chrono::system_clock::to_time_t(start);
	std::cout << "Time stamp: " << std::ctime(&startTime) << std::endl;

    /*! [Types] */
    using VectorDimI = LatticeCore<2>::VectorDimI;
    using IntScalarType = LatticeCore<2>::IntScalarType;
    /*! [Types] */

    /*! [Lattice] */
    const auto A(TextFileParser("bicrystal_2d.txt").readMatrix<double,2,2>("A",true));
    Lattice<2> lattice(A);
    /*! [Lattice] */

    /*! [Test] */
    const auto& coincidentLattices= lattice.generateCoincidentLattices(1e-3,10,15);
    /*! [Test] */

    /*! [SNF] */
    for (const auto& deformationGradient : coincidentLattices)
    {
        try
        {
            BiCrystal<2> bc(lattice,Lattice<2>(lattice.latticeBasis,deformationGradient),false);
            if (abs(bc.sigmaA) > 100000 || abs(bc.sigmaB) > 100000)
                continue;
            std::cout << "sigmaA = " << bc.sigmaA << "; sigmaB = " << bc.sigmaB << std::endl;

            Eigen::Matrix2d F= deformationGradient;
            std::cout << endl;
            std::cout << "Deformation gradient, F=" << std::endl;
            std::cout << F << std::endl;
            Eigen::Matrix2d C= F.transpose()*F;
            Eigen::Matrix2d E= (C-Eigen::Matrix2d::Identity())/2.0;

            std::cout << endl;
            std::cout << "Polar decomposition of F:" << std::endl;
            double s= sqrt(C(0,0)*C(1,1)-pow(C(0,1),2));
            double t= sqrt(C(0,0)+C(1,1)+2*s);
            Eigen::Matrix2d U= (C+s*Eigen::Matrix2d::Identity())/t;
            Eigen::Matrix2d R= F*U.inverse();
            std::cout << "R(" << atan2(R(1,0),R(0,0))*180/M_PI << ")= " << std::endl;
            std::cout << R << std::endl;
            std::cout << "U= " << std::endl;
            std::cout << U << std::endl;

            std::cout << endl;
            std::cout << "Elastic strain:" << std::endl;
            std::cout << E << std::endl;
            std::cout << "---------------------------------------------" << std::endl;
        }
        catch(std::runtime_error& e)
        {
            std::cout << e.what() << "Moving on ..." << std::endl;
        }
    }
    /*! [SNF] */

    auto end= std::chrono::system_clock::now();
    std::time_t endTime= std::chrono::system_clock::to_time_t(end);
    std::chrono::duration<double> elapsedSeconds(end-start);
    std::cout << "Elapsed time: " << elapsedSeconds.count() << " seconds" << std::endl;
    std::cout << "End of simulation" << std::endl;
    return 0;
}
