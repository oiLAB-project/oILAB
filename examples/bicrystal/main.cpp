/* This file is part of gbLAB.
 *
 * gbLAB is distributed without any warranty under the MIT License.
 */

#include <LatticeModule.h>
#include <TextFileParser.h>

using namespace gbLAB;

int main(int argc, char** argv)
{
    using IntScalarType= long long int;
    
    //const int case = argc > 1 ? std::atoi(argv[1]) : 2;
    const int caseId= 2;

    switch (caseId) {
        case 2:
        {
            const auto A(TextFileParser("bicrystal_2d.txt").readMatrix<double,2,2>("A",true));
            Lattice<2> L1(A);
//            Eigen::Matrix2d R= Eigen::Rotation2D<double>(38.04313507482*M_PI/180).matrix();
//            Eigen::Matrix2d R= Eigen::Rotation2D<double>(81.95686492519*M_PI/180).matrix();
//            Eigen::Matrix2d R= Eigen::Rotation2D<double>(38.05455500587*M_PI/180).matrix();
//            Eigen::Matrix2d R= Eigen::Rotation2D<double>(90.02872079759*M_PI/180).matrix();
//            Eigen::Matrix2d R= Eigen::Rotation2D<double>(43.89139748113*M_PI/180).matrix();
            Eigen::Matrix2d R= Eigen::Rotation2D<double>( 143.0029103206*M_PI/180).matrix();

            Lattice<2> L2(A,R);
            BiCrystal<2> bc(L1,L2);
            std::cout<<"sigma="<<bc.sigma<<std::endl;

            std::vector<LatticeVector<2>> boxVectors;
            Eigen::Matrix<IntScalarType,2,2> U;
            U = RLLL(bc.csl.latticeBasis, 0.75).unimodularMatrix();
            Eigen::Vector<IntScalarType,2> temp(U.col(0));
            boxVectors.push_back(LatticeVector<2>(temp,bc.csl));
            temp= U.col(1);
            boxVectors.push_back(LatticeVector<2>(temp,bc.csl));

            bc.box(boxVectors,0.6,"bc.txt",true);

            break;
        }
            
        case 3:
        {
            const auto A(TextFileParser("bicrystal_3d.txt").readMatrix<double,3,3>("A",true));
            const auto R1(TextFileParser("bicrystal_3d.txt").readMatrix<double,3,3>("R1",true));
            const auto R2(TextFileParser("bicrystal_3d.txt").readMatrix<double,3,3>("R2",true));
            const auto misAxis(TextFileParser("bicrystal_3d.txt").readMatrix<double,3,1>("misAxis",true));


            Lattice<3> L1(A,R1);
            Lattice<3> L2(A,R2);
            BiCrystal<3> bc(L1,L2);
            std::cout<<"sigma="<<bc.sigma<<std::endl;
            
            const auto a1(L1.latticeDirection(misAxis));

            std::cout<<"Integer coordinates of misorientation axis = "<<a1.latticeVector().transpose()<<std::endl;


            break;
        }
            
        case 4:
        {
            const auto A(TextFileParser("Grimmer1985.txt").readMatrix<double,3,3>("A",true));
            const auto Tn(TextFileParser("Grimmer1985.txt").readMatrix<long long int,3,3>("Tn",true));
            const auto Td(TextFileParser("Grimmer1985.txt").readScalar<long long int>("Td",true));
            RationalMatrix<3> T(Tn,Td);
            
            const Eigen::Matrix<double,3,3> R(A*T.asMatrix()*A.inverse());
            const double angle(acos(R.trace()/3.0));
            std::cout<<"rot angle="<<angle*180.0/M_PI<<std::endl;

            Eigen::Matrix<double,3,1> axis(Eigen::Matrix<double,3,1>::Zero());
            Eigen::EigenSolver< Eigen::Matrix<double,3,3> > es(R);

            for(int k=0;k<3;++k)
            {
                const auto& eVal(es.eigenvalues()(k));
                if(std::fabs(eVal.imag())<FLT_EPSILON && std::fabs(eVal.real()-1.0)<FLT_EPSILON)
                {
                    const auto eVec(es.eigenvectors().col(k));
                    for(int d=0;d<3;++d)
                    {
                        axis(d)=eVec(d).real();
                    }
                    break;
                }
            }
            std::cout<<"rot axis="<<axis.transpose()<<std::endl;

            
            Lattice<3> L1(A);
            const auto axisDir(L1.latticeDirection(axis));
            std::cout<<"axisDir="<<axisDir.latticeVector().transpose()<<std::endl;
            const auto axisRec(L1.reciprocalLatticeDirection(axis));
            std::cout<<"axisRec="<<axisRec.reciprocalLatticeVector().transpose()<<std::endl;

            
            Lattice<3> L2(A,R);
            BiCrystal<3> bc(L1,L2);
            std::cout<<"mu_R="<<bc.mu<<std::endl;
            std::cout<<"sigma_R="<<bc.sigma<<std::endl;
            
            std::cout<<"D=\n"<<bc.matrixD()<<std::endl;
            std::cout<<"M=\n"<<bc.M<<std::endl;
            std::cout<<"N=\n"<<bc.N<<std::endl;

            const Eigen::Matrix< long long,3,3> X(bc.matrixX());
            const Eigen::Matrix< long long,3,3> Y(bc.matrixY());
            const Eigen::Matrix< long long,3,3> adjX(X.adjoint());
            const Eigen::Matrix< long long,3,3> adjY(Y.adjoint());
            
            Lattice<3> L3(A,R.transpose());
            BiCrystal<3> bcT(L1,L3);
            std::cout<<"mu_RT="<<bcT.mu<<std::endl;
            std::cout<<"sigma_RT="<<bcT.sigma<<std::endl;

            std::cout<<"D=\n"<<bcT.matrixD()<<std::endl;
            std::cout<<"M=\n"<<bcT.M<<std::endl;
            std::cout<<"N=\n"<<bcT.N<<std::endl;

            break;
        }

        default:
        {
            std::cout<<"No input file for case="<< caseId <<std::endl;
            break;
        }
    }
    
    
    return 0;
}
