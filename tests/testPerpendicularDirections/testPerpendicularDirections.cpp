#include <IntegerLattice.h>
#include <iostream>
#include <random>

using namespace gbLAB;
int main()
{
    typedef long long int IntScalarType;

    const int dim=3;
    Eigen::Vector<IntScalarType,Eigen::Dynamic> in(dim);

    // generate random input direction
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::uniform_int_distribution<> distr(-10, 10); // define the range
    for(int n=0; n<dim; ++n)
    {
        in(n)= (IntScalarType) distr(gen);
    }

    in= in/IntegerMath<IntScalarType>::gcd(in);
    std::cout << "Input direction = " << in.transpose() << std::endl;

    // compute perpendicular output directions
    Eigen::Matrix<IntScalarType,Eigen::Dynamic,Eigen::Dynamic> out(dim-1,dim);
    out= IntegerLattice<dim>::perpendicularDirections(in);
    std::cout << "Perpendicular directions:" << std::endl;
    std::cout << out << std::endl;

    // Check if the directions are indeed perpendicular
    Eigen::Matrix<IntScalarType,Eigen::Dynamic,Eigen::Dynamic> mat(dim,dim);
    mat.block(1,0,dim-1,dim)= out;
    mat.row(0)= in;
    std::cout << "Dot products:" << std::endl;
    try {
        for (int i = 1; i < dim; i++) {
            std::cout << mat.row(0).dot(mat.row(i)) << std::endl;
            if (mat.row(0).dot(mat.row(i)) != 0)
                throw std::runtime_error("Output directions are not perpendicular to the given direction\n");
        }
    }
    catch (std::runtime_error& e){
        std::cout << e.what() << std::endl;
        return -1;
    }
    return 0;
}
