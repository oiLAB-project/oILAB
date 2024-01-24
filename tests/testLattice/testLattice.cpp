#include <LatticeModule.h>
#include <fstream>
#include<vector>

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
    std::cout << "Recover the reciprocal lattice vector from the reciprocal lattice direction s:" << std::endl;
    ReciprocalLatticeVector<3> t= s.reciprocalLatticeVector();
    std::cout << t <<  std::endl;
    /*! [Direction to vector] */

    /*! [Lattice direction along a cartesian vector] */
   LatticeDirection<3>  latticeDirectionAlong_s= L.latticeDirection(s.cartesian());
   std::cout << "Lattice direction along s = " << std::endl;
   std::cout << latticeDirectionAlong_s <<  std::endl;
    /*! [Lattice direction along a cartesian vector] */

    /*! [Stacking of a lattice plane] */
    std::cout << "Stacking of the lattice plane described by s = " << s.stacking() << std::endl;
    /*! [Stacking of a lattice plane] */

    /*! [Cross product] */
    ReciprocalLatticeDirection<3> uxv(u.cross(v));
    std::cout << "The cross product of the lattice vectors u and v is a lattice direction: " << uxv << std::endl;

    ReciprocalLatticeVector<3> q(L);
    q << 7,9,11;
    std::cout << "Reciprocal lattice vector q = " << std::endl;
    std::cout << q << std::endl;
    LatticeDirection<3> qxr(q.cross(r));
    std::cout << "The cross product of two reciprocal lattice vectors q and r is a lattice direction: " << qxr<< std::endl;
    /*! [Cross product] */


    /*! [Box3] */
    std::vector<LatticeVector<3>> boxVectors;
    boxVectors.push_back(LatticeVector<3>(u,L));
    boxVectors.push_back(LatticeVector<3>(v,L));
    boxVectors.push_back(LatticeVector<3>(latticeDirectionAlong_s.latticeVector(),L));
    std::cout << "Outputting a configuration of lattice points bounded by three box vectors: " << std::endl;
    for (const auto& vector : boxVectors)
        std::cout << vector.cartesian().transpose() << std::endl;
    L.box(boxVectors,"lattice.txt");
    /*! [Box3] */

    /*! [Box2] */
    Eigen::Matrix2d B;
    B << 1.0, 0.0,
         0.5, 1.0;
    Lattice<2> L2(B);
    LatticeVector<2> b1(L2); b1 << 1,13;
    LatticeVector<2> b2(L2); b2 << 2,7;
    std::cout << "Outputting a configuration of lattice points bounded by two box vectors: " << std::endl;
    std::cout << b1.cartesian().transpose() << std::endl;
    std::cout << b2.cartesian().transpose() << std::endl;
    L2.box(std::vector<LatticeVector<2>>{-1*b1,b2},"lattice2.txt");
    /*! [Box2] */

    return 0;
}
