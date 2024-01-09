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
    //const auto& coincidentLattices= lattice.generateCoincidentLattices(1e-3,10,15);
    const auto& coincidentLattices= lattice.generateCoincidentLattices(0,300,30000);
    //const auto& coincidentLattices= lattice.generateCoincidentLattices(1e-2,30,15);
    /*! [Test] */

    /*! [SNF] */
    for (const auto& deformationGradient : coincidentLattices)
    {
        try
        {
            BiCrystal<2> bc(lattice,Lattice<2>(lattice.latticeBasis,deformationGradient),false);
            if (abs(bc.sigmaA) > 60000 || abs(bc.sigmaB) > 60000)
                continue;
            std::cout << std::endl;
            std::cout << "--------------------------------- SNF -------------------------------------" << std::endl;
            std::cout << "sigmaA = " << bc.sigmaA << "; sigmaB = " << bc.sigmaB << std::endl;
            std::cout << "Ap= " << std::endl;
            std::cout << std::setprecision(12) << bc.Ap.latticeBasis << std::endl;
            std::cout << "Bp= " << std::endl;
            std::cout << std::setprecision(12) << bc.Bp.latticeBasis << std::endl;
            std::cout << "CSL= " << std::endl;
            std::cout << std::setprecision(12) << bc.csl.latticeBasis << std::endl;
            std::cout << "DCSL= " << std::endl;
            std::cout << std::setprecision(12) << bc.dscl.latticeBasis << std::endl;
            std::cout << endl;
            /*! [SNF] */

            /*! [Heterodeform] */
            std::cout << "------Hetero-deformation -------" << std::endl;
            Eigen::Matrix2d F= deformationGradient;
            std::cout << "Deformation gradient, F=" << std::endl;
            std::cout << std::setprecision(12) << F << std::endl;
            Eigen::Matrix2d C= F.transpose()*F;
            Eigen::Matrix2d E= (C-Eigen::Matrix2d::Identity())/2.0;

            std::cout << endl;
            std::cout << "Polar decomposition of F:" << std::endl;
            double s= sqrt(C(0,0)*C(1,1)-pow(C(0,1),2));
            double t= sqrt(C(0,0)+C(1,1)+2*s);
            Eigen::Matrix2d U= (C+s*Eigen::Matrix2d::Identity())/t;
            Eigen::Matrix2d R= F*U.inverse();
            std::cout << "R(" << std::setprecision(12) << atan2(R(1,0),R(0,0))*180/M_PI << ")= " << std::endl;
            std::cout << R << std::endl;
            std::cout << "U= " << std::endl;
            std::cout << std::setprecision(12) << U << std::endl;

            std::cout << endl;
            std::cout << "Elastic strain:" << std::endl;
            std::cout << std::setprecision(12) << E << std::endl;
            std::cout << endl;
            /*! [Heterodeform] */

            /*! [Invariance] */
            std::cout << "------Invariance of the moire -------" << std::endl;
            auto reducedDsclBasis= RLLL(bc.dscl.latticeBasis,0.75);
            auto U_Dscl= reducedDsclBasis.unimodularMatrix();

            std::cout << "Reduced DSCL basis vectors:" << std::endl;
            std::cout << "d1 = ";
            std::cout << std::setprecision(20) << reducedDsclBasis.reducedBasis().col(0).transpose() << std::endl;
            std::cout << "Integer coordinates of d1:";
            LatticeVector<2> d1(bc.dscl);
            d1 << U_Dscl.col(0).template cast<IntScalarType>();
            std::cout << std::setprecision(20) << d1.transpose() << std::endl;
            std::cout << std::endl;

            LatticeVector<2> d2(bc.dscl);
            std::cout << "d2 = ";
            std::cout << std::setprecision(20) << reducedDsclBasis.reducedBasis().col(1).transpose() << std::endl;
            std::cout << "Integer coordinates of d2:";
            d2 << U_Dscl.col(1).template cast<IntScalarType>();
            std::cout << std::setprecision(20) << d2.transpose() << std::endl;
            std::cout << std::endl;

            // shift vectors corresponding to d1 and d2
            LatticeVector<2> s1(bc.dscl), s2(bc.dscl);
            s1 << bc.LambdaA * d1;
            s2 << bc.LambdaA * d2;
            Lattice<2> reducedCsl(RLLL(bc.csl.latticeBasis,0.75).reducedBasis());
            std::cout << "Reduced shift vectors: " << std::endl;

            Eigen::Vector2d s1_coordinates_in_reduced_csl= reducedCsl.latticeBasis.inverse()*s1.cartesian();
            Eigen::Vector2d s2_coordinates_in_reduced_csl= reducedCsl.latticeBasis.inverse()*s2.cartesian();
            Eigen::Vector2d s1_coordinates_modulo= s1_coordinates_in_reduced_csl.array()-s1_coordinates_in_reduced_csl.array().round();
            Eigen::Vector2d s2_coordinates_modulo= s2_coordinates_in_reduced_csl.array()-s2_coordinates_in_reduced_csl.array().round();
            std::cout << "s1 = ";
            std::cout << std::setprecision(20) << (reducedCsl.latticeBasis * s1_coordinates_modulo).transpose() << std::endl;
            std::cout << "s2 = ";
            std::cout << std::setprecision(20) << (reducedCsl.latticeBasis * s2_coordinates_modulo).transpose() << std::endl;
            std::cout << "---------------------------------------------------------------------------" << std::endl;
            /*! [Invariance] */

        }
        catch(std::runtime_error& e)
        {
            std::cout << e.what() << "Moving on ..." << std::endl;
        }
    }

    auto end= std::chrono::system_clock::now();
    std::time_t endTime= std::chrono::system_clock::to_time_t(end);
    std::chrono::duration<double> elapsedSeconds(end-start);
    std::cout << "Elapsed time: " << elapsedSeconds.count() << " seconds" << std::endl;
    std::cout << "End of simulation" << std::endl;
    return 0;
}
