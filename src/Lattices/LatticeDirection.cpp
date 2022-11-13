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
    LatticeDirection<dim>::LatticeDirection(const VectorDimI& v,
                     const Lattice<dim>& lat) :
    /* base init */ LatticeVector<dim>(((v.squaredNorm()==0)? v : (v/IntegerMath<IntScalarType>::gcd(v)).eval()),lat)
    {
    }

    template<int dim>
    basic_ostream<char>& operator<<(basic_ostream<char>& s, const LatticeDirection<dim>& m)
    {
        return s<<m.latticeVector().transpose();
    }

    template struct LatticeDirection<1>;
    template basic_ostream<char>& operator<<(basic_ostream<char>& s, const LatticeDirection<1>& m);
    template struct LatticeDirection<2>;
    template basic_ostream<char>& operator<<(basic_ostream<char>& s, const LatticeDirection<2>& m);
    template struct LatticeDirection<3>;
    template basic_ostream<char>& operator<<(basic_ostream<char>& s, const LatticeDirection<3>& m);
    template struct LatticeDirection<4>;
    template basic_ostream<char>& operator<<(basic_ostream<char>& s, const LatticeDirection<4>& m);
    template struct LatticeDirection<5>;
    template basic_ostream<char>& operator<<(basic_ostream<char>& s, const LatticeDirection<5>& m);

} // end namespace
#endif
