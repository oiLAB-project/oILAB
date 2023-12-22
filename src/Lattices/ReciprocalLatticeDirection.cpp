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

    template<int dim>
    basic_ostream<char>& operator<<(basic_ostream<char>& s, const ReciprocalLatticeDirection<dim>& m)
    {
        return s<<m.reciprocalLatticeVector().transpose();
    }

    template <int dim>
    double ReciprocalLatticeDirection<dim>::planeSpacing() const
    {
        return 1.0 / cartesian().norm();
    }

    template <int dim>
    int ReciprocalLatticeDirection<dim>::stacking() const
    {
        VectorDimD doubleCoordinates= lattice.reciprocalBasis.transpose()*cartesian();
        VectorDimI integerCoordinates= LatticeCore<dim>::rationalApproximation(doubleCoordinates);
        LatticeDirection<dim> directionAlongReciprocalDirection(integerCoordinates,lattice);
        return dot(directionAlongReciprocalDirection);
return 1;
    }

    template struct ReciprocalLatticeDirection<1>;
    template basic_ostream<char>& operator<<(basic_ostream<char>& s, const ReciprocalLatticeDirection<1>& m);
    template struct ReciprocalLatticeDirection<2>;
    template basic_ostream<char>& operator<<(basic_ostream<char>& s, const ReciprocalLatticeDirection<2>& m);
    template struct ReciprocalLatticeDirection<3>;
    template basic_ostream<char>& operator<<(basic_ostream<char>& s, const ReciprocalLatticeDirection<3>& m);
    template struct ReciprocalLatticeDirection<4>;
    template basic_ostream<char>& operator<<(basic_ostream<char>& s, const ReciprocalLatticeDirection<4>& m);
    template struct ReciprocalLatticeDirection<5>;
    template basic_ostream<char>& operator<<(basic_ostream<char>& s, const ReciprocalLatticeDirection<5>& m);

} // end namespace
#endif
