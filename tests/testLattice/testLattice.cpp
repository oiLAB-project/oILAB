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
    std::cout << "Recover lattice vector from its cartesian components:";
    std::cout << L.latticeVector(u.cartesian()) << std::endl;
    /*! [Integer coordinates] */
    /*! [Cartesian coordinates] */
    std::cout << "Create a lattice vector using cartesian coordinates: 1.5,2,3.5" << std::endl;
    Eigen::Vector3d vCoords; vCoords << 1.5,2,3.5;
    LatticeVector<3> v(vCoords, L);
    std::cout << "Integer coordinates of v:" << std::endl;
    std::cout << v << std::endl;
    /*! [Cartesian coordinates] */
    /*! [Cartesian coordinates fail] */
    try{
        // This fails
        Eigen::Vector3d wCoords; wCoords << 1.8,2,3.5;
        LatticeVector<3> w(wCoords, L);
    }
    catch(std::runtime_error& e)
    {
        std::cout << e.what() << std::endl;
    }
    /*! [Cartesian coordinates fail] */

    /*! [Lattice vector algebra] */
    v << 1,1,2;
    LatticeVector<3> w= 2*u + 6*v;
    std::cout << "Integer coordinates of w:" << std::endl;
    std::cout << w << std::endl;
    std::cout << "Cartesian coordinates of w:" << std::endl;
    std::cout << w.cartesian() << std::endl;
    /*! [Lattice vector algebra] */

    /*! [Lattice direction] */
    std::cout << "Lattice direction along w:" << std::endl;
    LatticeDirection<3> wDir(w);
    std::cout << wDir << std::endl;
    /*! [Lattice direction] */


    /*! [Reciprocal lattice vector] */
    std::cout << "Creating a reciprocal lattice vector r:" << std::endl;
    ReciprocalLatticeVector<3> r(L);
    r << 3,6,9;
    std::cout << r << std::endl;
    std::cout << "Cartesian coordinates of the reciprocal lattice vector r:" << std::endl;
    std::cout << r.cartesian() << std::endl;
    /*! [Reciprocal lattice vector] */

    /*! [Reciprocal lattice direction] */
    std::cout << "Reciprocal lattice direction along r:" << std::endl;
    ReciprocalLatticeDirection<3> s(r);
    std::cout << s <<  std::endl;
    /*! [Reciprocal lattice direction] */

    /*! [Direction to vector] */
    std::cout << "Get reciprocal lattice vector from reciprocal lattice direction:" << std::endl;
    ReciprocalLatticeVector<3> t= s.reciprocalLatticeVector();
    std::cout << t <<  std::endl;
    /*! [Direction to vector] */

    /*! [Direction along a reciprocal direction] */
    std::cout << "Get a lattice direction along the reciprocal lattice direction " << s << std::endl;
    Eigen::Vector3d doubleCoordinates= L.reciprocalBasis.transpose()*s.cartesian();
    auto integerCoordinates= LatticeCore<3>::rationalApproximation(doubleCoordinates);
    LatticeDirection<3> directionAlongReciprocalDirection(integerCoordinates,L);
    std::cout << directionAlongReciprocalDirection <<  std::endl;
    std::cout << "Stacking for the above plane = " << directionAlongReciprocalDirection.dot(s) << std::endl;
    /*! [Direction along a reciprocal direction] */

    /*! [Cross product] */
   ReciprocalLatticeDirection<3> uxv(u.cross(v));
   LatticeDirection<3> rxs(r.cross(s.reciprocalLatticeVector()));
    /*! [Cross product] */
    return 0;
}
