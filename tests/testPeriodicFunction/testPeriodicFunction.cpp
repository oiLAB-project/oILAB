#include <Eigen/Core>
#include <LatticeModule.h>
#include <PeriodicFunction.h>
#include <fstream>

using namespace gbLAB;

#define MY_ERROR(message)                                                \
  {                                                                      \
    std::cout << "* Error : \"" << message << "\" : " << __LINE__ << ":" \
              << __FILE__ << std::endl;                                  \
    exit(1);                                                             \
  }

int main()
{
    // Read the kernel function data
    std::map<double,dcomplex> values;
    double tmp1; dcomplex tmp2;
    std::ifstream file3D("meyer.txt");
    if(!file3D) MY_ERROR(std::string("ERROR: meyer.txt could not be opened for reading!"));
    for (int i=0; i<16001;i++) {
        file3D >> tmp1 >> tmp2;
        values[tmp1]= tmp2;
    }

    // create the kernel function
    PiecewisePolynomial<dcomplex,2> kernelFunction(values);

    // Lattice system
    Eigen::array<Eigen::Index,2> n{32,32};
    Eigen::Matrix2d A;
    A << -1.23,          2.46,
         -2.13042249331, 0.00;
    Lattice<2> L(A);
    Eigen::Matrix<double,2,2> basisAtoms;
    basisAtoms.col(0).setZero();
    basisAtoms.col(1)= A.col(0)/3 + 2*A.col(1)/3;
    MultiLattice<2> M(A,basisAtoms);

    // Generate potential V
    PeriodicFunction<dcomplex ,2> V2D= PeriodicFunction<dcomplex,2>::kernelConvolution(n,M,kernelFunction);
    std::cout << V2D.values << std::endl;

    Eigen::Vector2d center= A.colwise().sum()/2.0;
    double potential= 0;
    for (int i=-8; i<8; i++)
    {
        for (int j=-8; j<8; j++)
        {
            for (const auto& atom : basisAtoms.colwise()) {
                Eigen::Vector2d position = i * A.col(0) + j * A.col(1) + atom;
                potential+= kernelFunction(center-position).real();
            }
        }
    }
    std::cout << potential << std::endl;
    std::cout << V2D.values(16,16)<< std::endl;


    return 0;
}

