#include <TextFileParser.h>
#include <LatticeModule.h>
#include <random>

using namespace gbLAB;
int main()
{
    using IntScalarType= long long int;

    /*! [Random] */
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::uniform_int_distribution<> distr(-10, 10); // define the range
    /*! [Random] */

    /*! [Lattice] */
    const int dim= 5;
    Eigen::Matrix<double,dim,dim> A;
    A.setIdentity();
    std::cout << "Testing in dimension = " << dim << std::endl;
    Lattice<dim> lat(A);
    /*! [Lattice] */

    /*! [Directions] */
    Eigen::Vector<IntScalarType,dim> millerIndices;
    for (auto& element : millerIndices)
        element= distr(gen);
    ReciprocalLatticeDirection<dim> rDir(ReciprocalLatticeVector<dim>(millerIndices,lat));

    Eigen::Vector<IntScalarType,dim> latticeCoordinates;
    for (auto& element : latticeCoordinates)
        element= distr(gen);
    LatticeDirection<dim> lDir(LatticeVector<dim>(latticeCoordinates,lat));
    /*! [Directions] */

    std::cout << "Input Miller index:" << rDir.transpose() << std::endl;
    /*! [Basis1] */
    auto directions=  lat.planeParallelLatticeBasis(rDir);
    std::cout << "Plane parallel lattice basis: " << std::endl;
    for (auto it= directions.begin(); it!=directions.end(); ++it)
    {
        std::cout << it->transpose() << std::endl;
    }
    /*! [Basis1] */
    /*! [Test1] */
    try {
        for (auto it = directions.begin(); it != directions.end(); ++it) {
            if (it == directions.begin())
            {
                if (rDir.dot(*it) != 1)
                    throw std::runtime_error("The dot product of the input reciprocal vector and the first basis vector is not +1 or -1\n");
            }
            else
            {
                if (rDir.dot(*it) != 0)
                    throw std::runtime_error(
                            "One of the plane-parallel basis vector is not perpendicular to the input reciprocal lattice direction\n");
            }
        }
    }
    /*! [Test1] */
    catch (std::runtime_error& e){
        std::cout << e.what() << std::endl;
        return -1;
    }

    std::cout << "----------------------------------------------------" << std::endl;

    // test the member function directionOrthogonalReciprocalLatticeBasis function
    std::cout << "Input lattice direction :" << lDir.transpose() << std::endl;
    /*! [Basis2] */
    auto reciprocalDirections=  lat.directionOrthogonalReciprocalLatticeBasis(lDir);
    std::cout << "Direction orthogonal reciprocal basis: " << std::endl;
    for (auto it= reciprocalDirections.begin(); it!=reciprocalDirections.end(); ++it)
    {
        std::cout << it->transpose() << std::endl;
    }
    /*! [Basis2] */
    /*! [Test2] */
    try {
        for (auto it = reciprocalDirections.begin(); it != reciprocalDirections.end(); ++it) {
            if (it == reciprocalDirections.begin())
            {
                if (lDir.dot(*it) != 1)
                    throw std::runtime_error("The dot product of the input lattice direction and the first reciprocal basis vector is not +1 or -1\n");
            }
            else
            {
                if (lDir.dot(*it) != 0)
                    throw std::runtime_error(
                            "One of the direction-orthogonal reciprocal basis vector is not perpendicular to the input lattice direction\n");
            }
        }
    }
    /*! [Test2] */
    catch (std::runtime_error& e){
        std::cout << e.what() << std::endl;
        return -1;
    }
    return 0;
}
