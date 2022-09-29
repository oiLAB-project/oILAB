#include <TextFileParser.h>
#include <LatticeModule.h>
#include <random>

using namespace gbLAB;
int main()
{
    using IntScalarType= long long int;
    const auto A(TextFileParser("bicrystal_2d.txt").readMatrix<double,2,2>("A",true));

    // setup the random number distribution for generating random input reciprocal directions
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::uniform_int_distribution<> distr(-10, 10); // define the range

    // Testing for two dimensions
    std::cout << "Testing for two dimensions:" << std::endl;
    Lattice<2> L2(A);
    Eigen::Vector<IntScalarType ,Eigen::Dynamic> d2(2);
    d2 << distr(gen),distr(gen);
    ReciprocalLatticeDirection<2> ld2(ReciprocalLatticeVector<2>(d2,L2));
    std::cout << "Input Miller index:" << ld2.transpose() << std::endl;

    auto directions=  L2.planeParallelLatticeDirections(ld2);
    std::cout << "Plane parallel directions:" << std::endl;
    for (auto it= directions.begin(); it!=directions.end(); ++it) {
        std::cout << it->transpose() << std::endl;
        if (ld2.dot(*it) != 0) throw std::runtime_error("Output directions are not perpendicular to the given reciprocal vector\n");
    }



    // Testing for three dimensions
    std::cout << "Testing for three dimensions:" << std::endl;
    const auto B(TextFileParser("bicrystal_3d.txt").readMatrix<double,3,3>("B",true));
    const auto R(TextFileParser("bicrystal_3d.txt").readMatrix<double,3,3>("R",true));
    Lattice<3> L3(B,R);
    Eigen::Vector<IntScalarType ,Eigen::Dynamic> d3(3);
    d3 << distr(gen),distr(gen),distr(gen);
    ReciprocalLatticeDirection<3> ld3(ReciprocalLatticeVector<3>(d3,L3));
    std::cout << "Input Miller index:" << ld3.transpose() << std::endl;

    std::cout << "Plane parallel directions:" << std::endl;
    auto directions3=  L3.planeParallelLatticeDirections(ld3);
    for (auto it= directions3.begin(); it!=directions3.end(); ++it) {
        std::cout << it->transpose() << std::endl;
        if (ld3.dot(*it) != 0)
            throw std::runtime_error("Output directions are not perpendicular to the given reciprocal vector\n");
    }

    return 0;
}
