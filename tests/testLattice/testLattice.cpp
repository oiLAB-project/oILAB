#include <LatticeModule.h>

using namespace gbLAB;
int main()
{
    /*! [Lattice] */
    Eigen::Matrix3d A;
    A << 0.5, 0.5, 0,
         0.5, 0.0, 0.5,
         0.0, 0.5, 0.5;
    Lattice<3> L(A);

    std::cout << "Lattice basis:" << std::endl;
    std::cout << L.latticeBasis << std::endl;
    std::cout << "Reciprocal lattice basis:" << std::endl;
    std::cout << L.reciprocalBasis << std::endl;
    /*! [Lattice] */

    /*! [Lattice vector] */
    std::cout << "Creating a zero lattice vector u:" << std::endl;
    LatticeVector<3> u(L);
    std::cout << u << std::endl;
    /*! [Lattice vector] */
    /*! [Integer coordinates] */
    std::cout << "Modifying u:" << std::endl;
    u << 2,1,3;
    std::cout << u << std::endl;
    std::cout << "Cartesian coordinates of u:" << std::endl;
    std::cout << u.cartesian() << std::endl;
    /*! [Integer coordinates] */
    /*! [Cartesian coordinates] */
    std::cout << "Create a lattice vector using cartesian coordinates: 1.5,2,3.5" << std::endl;
    Eigen::Vector3d vCoords; vCoords << 1.5,2,3.5;
    LatticeVector<3> v(vCoords,L);
    std::cout << "Integer coordinates of v:" << std::endl;
    std::cout << v << std::endl;
    /*! [Cartesian coordinates] */

    /*! [Lattice vector algebra] */
    v << 1,1,2;
    LatticeVector<3> w= 2*u + 3*v;
    std::cout << "Cartesian coordinates of the lattice vector w:" << std::endl;
    std::cout << w.cartesian() << std::endl;
    /*! [Lattice vector algebra] */

    /*! [Reciprocal lattice vector ] */
    // Reciprocal lattice vector
    std::cout << "Creating a zero reciprocal lattice vector u:" << std::endl;
    ReciprocalLatticeVector<3> r(L);
    std::cout << r << std::endl;
    std::cout << "Cartesian coordinates of the reciprocal lattice vector u:" << std::endl;
    std::cout << r.cartesian() << std::endl;
    /*! [Reciprocal lattice vector ] */
    return 0;
}
