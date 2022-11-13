/* This file is part of gbLAB.
 *
 * gbLAB is distributed without any warranty under the MIT License.
 */


#ifndef gbLAB_ReciprocalLatticeDirection_cpp_
#define gbLAB_ReciprocalLatticeDirection_cpp_

#include <LatticeModule.h>

namespace gbLAB
{
    template <int dim>
    ReciprocalLatticeDirection<dim>::ReciprocalLatticeDirection(const ReciprocalLatticeVector<dim>& v) :
    /* init */ ReciprocalLatticeVector<dim>(((v.squaredNorm()==0)? v : (v/IntegerMath<IntScalarType>::gcd(v)).eval()),v.lattice)
    {
    }


    template struct ReciprocalLatticeDirection<1>;
    template struct ReciprocalLatticeDirection<2>;
    template struct ReciprocalLatticeDirection<3>;
    template struct ReciprocalLatticeDirection<4>;
    template struct ReciprocalLatticeDirection<5>;

} // end namespace
#endif
