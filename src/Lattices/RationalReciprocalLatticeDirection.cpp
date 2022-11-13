/* This file is part of gbLAB.
 *
 * gbLAB is distributed without any warranty under the MIT License.
 */


#ifndef gbLAB_RationalReciprocalLatticeDirection_cpp_
#define gbLAB_RationalReciprocalLatticeDirection_cpp_

#include <LatticeModule.h>
#include <RationalReciprocalLatticeDirection.h>

namespace gbLAB
{
    /**********************************************************************/
    template <int dim>
    RationalReciprocalLatticeDirection<dim>::RationalReciprocalLatticeDirection(const Rational<IntScalarType> &_rat, const ReciprocalLatticeDirection<dim> &_dir) : /* init */ rat(_rat)
                                                                                        /* init */,
                                                                                        dir(_dir)
    {
    }

    /**********************************************************************/
    template <int dim>
    RationalReciprocalLatticeDirection<dim>::RationalReciprocalLatticeDirection(const Rational<IntScalarType> &_rat, const ReciprocalLatticeVector<dim> &v) :
    /* init */ RationalReciprocalLatticeDirection(_rat, ReciprocalLatticeDirection<dim>(v))
    {
    }

    /**********************************************************************/
    template <int dim>
    RationalReciprocalLatticeDirection<dim>::RationalReciprocalLatticeDirection(const ReciprocalLatticeVector<dim> &v) :
    /* init */ rat(Rational<IntScalarType>(IntegerMath<IntScalarType>::gcd(v), 1)),
    /* init */ dir(v)
    {
    }

    /**********************************************************************/

    template <int dim>
    typename RationalReciprocalLatticeDirection<dim>::VectorDimD RationalReciprocalLatticeDirection<dim>::cartesian() const
    {
        return dir.cartesian() * rat.asDouble();
    }

    /**********************************************************************/
    template <int dim>
    Rational<typename RationalReciprocalLatticeDirection<dim>::IntScalarType> RationalReciprocalLatticeDirection<dim>::dot(const LatticeVector<dim> &other) const
    {
        return rat * other.dot(dir);
    }

    /**********************************************************************/
    template <int dim>
    RationalReciprocalLatticeDirection<dim> RationalReciprocalLatticeDirection<dim>::operator*(const IntScalarType &scalar) const
    {
        return RationalReciprocalLatticeDirection<dim>(rat * scalar, dir);
    }

    /**********************************************************************/
    template <int dim>
    RationalReciprocalLatticeDirection<dim> RationalReciprocalLatticeDirection<dim>::operator/(const IntScalarType &scalar) const
    {
        return RationalReciprocalLatticeDirection<dim>(rat / scalar, dir);
    }

    /**********************************************************************/
    template <int dim>
    RationalReciprocalLatticeDirection<dim> RationalReciprocalLatticeDirection<dim>::operator+(const RationalReciprocalLatticeDirection<dim> &other) const
    {
        assert(&dir.lattice == &other.dir.lattice && "Rational Lattice Vector Type belong to different Lattices.");
        //const VectorDimI temp(rat.n * other.rat.d * dir + other.rat.n * rat.d * other.dir);
        VectorDimI temp(rat.n * other.rat.d * dir + other.rat.n * rat.d * other.dir);
        const IntScalarType gcd(IntegerMath<IntScalarType>::gcd(temp));
        temp= temp/gcd;
        //const ReciprocalLatticeVector<dim> v(temp / gcd, dir.lattice);
        const ReciprocalLatticeVector<dim> v(temp, dir.lattice);
        return RationalReciprocalLatticeDirection<dim>(Rational(gcd, rat.d * other.rat.d), ReciprocalLatticeDirection<dim>(v));
    }

    /**********************************************************************/
    template <int dim>
    RationalReciprocalLatticeDirection<dim> RationalReciprocalLatticeDirection<dim>::operator-(const RationalReciprocalLatticeDirection<dim> &other) const
    {
        assert(&dir.lattice == &other.dir.lattice && "ReciprocalLatticeVectorType belong to different Lattices.");
        //const VectorDimI temp(rat.n * other.rat.d * dir - other.rat.n * rat.d * other.dir);
        VectorDimI temp(rat.n * other.rat.d * dir - other.rat.n * rat.d * other.dir);
        const IntScalarType gcd(IntegerMath<IntScalarType>::gcd(temp));
        temp= temp/gcd;
        //const ReciprocalLatticeVector<dim> v(temp / gcd, dir.lattice);
        const ReciprocalLatticeVector<dim> v(temp, dir.lattice);
        return RationalReciprocalLatticeDirection<dim>(Rational(gcd, rat.d * other.rat.d), ReciprocalLatticeDirection<dim>(v));
    }

    /**********************************************************************/
    template <int dim>
    RationalReciprocalLatticeDirection<dim> RationalReciprocalLatticeDirection<dim>::operator+(const ReciprocalLatticeVector<dim> &other) const
    {
        assert(&dir.lattice == &other.lattice && "ReciprocalLatticeVectorType belong to different Lattices.");
        const IntScalarType gcd(IntegerMath<IntScalarType>::gcd(other));
        return this->operator+(RationalReciprocalLatticeDirection<dim>(Rational<IntScalarType>(gcd, 1), ReciprocalLatticeDirection<dim>(other)));
    }

    /**********************************************************************/
    template <int dim>
    RationalReciprocalLatticeDirection<dim> RationalReciprocalLatticeDirection<dim>::operator-(const ReciprocalLatticeVector<dim> &other) const
    {
        assert(&dir.lattice == &other.lattice && "ReciprocalLatticeVectorType belong to different Lattices.");
        const IntScalarType gcd(IntegerMath<IntScalarType>::gcd(other));
        return this->operator-(RationalReciprocalLatticeDirection<dim>(Rational<IntScalarType>(gcd, 1), ReciprocalLatticeDirection<dim>(other)));
    }

    /**********************************************************************/
    template <int dim>
    double RationalReciprocalLatticeDirection<dim>::squaredNorm() const
    {
        return dir.squaredNorm() * std::pow(rat.asDouble(), 2);
    }
    
    template<int dim>
    RationalReciprocalLatticeDirection<dim> operator*(const typename RationalReciprocalLatticeDirection<dim>::IntScalarType& scalar, const RationalReciprocalLatticeDirection<dim>& L)
    {
        return L*scalar;
    }
    
    template struct RationalReciprocalLatticeDirection<1>;
    template RationalReciprocalLatticeDirection<1> operator*(const typename RationalReciprocalLatticeDirection<1>::IntScalarType& scalar, const RationalReciprocalLatticeDirection<1>& L);
    template struct RationalReciprocalLatticeDirection<2>;
    template RationalReciprocalLatticeDirection<2> operator*(const typename RationalReciprocalLatticeDirection<2>::IntScalarType& scalar, const RationalReciprocalLatticeDirection<2>& L);
    template struct RationalReciprocalLatticeDirection<3>;
    template RationalReciprocalLatticeDirection<3> operator*(const typename RationalReciprocalLatticeDirection<3>::IntScalarType& scalar, const RationalReciprocalLatticeDirection<3>& L);
    template struct RationalReciprocalLatticeDirection<4>;
    template RationalReciprocalLatticeDirection<4> operator*(const typename RationalReciprocalLatticeDirection<4>::IntScalarType& scalar, const RationalReciprocalLatticeDirection<4>& L);
    template struct RationalReciprocalLatticeDirection<5>;
    template RationalReciprocalLatticeDirection<5> operator*(const typename RationalReciprocalLatticeDirection<5>::IntScalarType& scalar, const RationalReciprocalLatticeDirection<5>& L);

} // end namespace
#endif
