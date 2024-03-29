/* This file is part of gbLAB.
 *
 * gbLAB is distributed without any warranty under the MIT License.
 */


#ifndef gbLAB_LatticeMath_
#define gbLAB_LatticeMath_

namespace gbLAB
{

    template <int dim>
    class Lattice;

    template <int dim>
    class LatticeVector;

    template <int dim>
    class ReciprocalLatticeVector;

    template <int dim>
    struct LatticeDirection;
    
    template <int dim>
    struct ReciprocalLatticeDirection;

    template <int dim>
    struct RationalLatticeDirection;

    template <int dim>
    struct RationalReciprocalLatticeDirection;

    template <int dim>
    class BiCrystal ;

    template <int dim>
    class Gb;

//    struct LatticePlane;
//    struct LatticeLine;

} // end namespace

//#include <LatticeGCD.h>
//#include <LatticeBase.h>
#include <LatticeCore.h>
#include <Lattice.h>
#include <LatticeVector.h>
#include <ReciprocalLatticeVector.h>
#include <LatticeDirection.h>
#include <ReciprocalLatticeDirection.h>
#include <RationalLatticeDirection.h>
#include <RationalReciprocalLatticeDirection.h>
#include <RationalMatrix.h>
#include <SmithDecomposition.h>
#include <BiCrystal.h>
#include <Gb.h>

#endif
