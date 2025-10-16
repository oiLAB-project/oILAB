/* This file is part of gbLAB.
 *
 * gbLAB is distributed without any warranty under the MIT License.
 */


#ifndef gbLAB_LatticeMath_
#define gbLAB_LatticeMath_

namespace oILAB {

template <int dim> class Lattice;

template <int dim> class LatticeVector;

template <int dim> class ReciprocalLatticeVector;

template <int dim> struct LatticeDirection;

template <int dim> struct ReciprocalLatticeDirection;

template <int dim> struct RationalLatticeDirection;

template <int dim> struct RationalReciprocalLatticeDirection;

template <int dim> class BiCrystal;

template <int dim> class Gb;

} // namespace oILAB

#include "LatticeCore.h"
#include "../Math/RationalMatrix.h"
#include "BiCrystal.h"
#include "Gb.h"
#include "Lattice.h"
#include "LatticeDirection.h"
#include "LatticeVector.h"
#include "RationalLatticeDirection.h"
#include "RationalReciprocalLatticeDirection.h"
#include "ReciprocalLatticeDirection.h"
#include "ReciprocalLatticeVector.h"
#include "../Math/SmithDecomposition.h"

#endif
