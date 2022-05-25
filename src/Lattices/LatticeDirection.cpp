/* This file is part of gbLAB.
 *
 * gbLAB is distributed without any warranty under the MIT License.
 */


#ifndef gbLAB_LatticeDirection_cpp_
#define gbLAB_LatticeDirection_cpp_

#include <LatticeModule.h>
namespace gbLAB
{

    
    template <int dim>
    LatticeDirection<dim>::LatticeDirection(const LatticeVector<dim>& v) :
    /* base init */ LatticeVector<dim>(((v.squaredNorm()==0)? v : (v/IntegerMath<IntScalarType>::gcd(v)).eval()),v.lattice)
    {
    }
        

    
    template struct LatticeDirection<1>;
    template struct LatticeDirection<2>;
    template struct LatticeDirection<3>;

} // end namespace
#endif
