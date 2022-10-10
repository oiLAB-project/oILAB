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
        
    template <int dim>
    LatticeDirection<dim>::LatticeDirection(const Eigen::Matrix<IntScalarType,dim,1>& v,
                     const Lattice<dim>& lat) :
    /* base init */ LatticeVector<dim>(((v.squaredNorm()==0)? v : (v/IntegerMath<IntScalarType>::gcd(v)).eval()),lat)
    {
    }

    template struct LatticeDirection<1>;
    template struct LatticeDirection<2>;
    template struct LatticeDirection<3>;
    template struct LatticeDirection<4>;
    template struct LatticeDirection<5>;

} // end namespace
#endif
