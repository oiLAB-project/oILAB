/* This file is part of gbLAB.
 *
 * gbLAB is distributed without any warranty under the MIT License.
 */


#ifndef gbLAB_RationalReciprocalLatticeDirection_h_
#define gbLAB_RationalReciprocalLatticeDirection_h_

#include <LatticeVector.h>
#include <ReciprocalLatticeVector.h>
#include <Rational.h>
#include <BestRationalApproximation.h>

namespace gbLAB
{
    template <int dim>
    struct RationalReciprocalLatticeDirection
    {
        typedef typename LatticeCore<dim>::IntScalarType IntScalarType;
        typedef typename LatticeCore<dim>::VectorDimD VectorDimD;
        typedef typename LatticeCore<dim>::MatrixDimD MatrixDimD;
        typedef typename LatticeCore<dim>::VectorDimI VectorDimI;
        typedef typename LatticeCore<dim>::MatrixDimI MatrixDimI;

        const Rational<IntScalarType> rat;
        const ReciprocalLatticeDirection<dim> dir;
        
    public:
        
        
        RationalReciprocalLatticeDirection(const Rational<IntScalarType>& _rat, const ReciprocalLatticeDirection<dim>& _dir) ;
        RationalReciprocalLatticeDirection(const Rational<IntScalarType>& _rat, const ReciprocalLatticeVector<dim>& v) ;
        RationalReciprocalLatticeDirection(const ReciprocalLatticeVector<dim>& v) ;
        RationalReciprocalLatticeDirection(const RationalReciprocalLatticeDirection<dim>& other) = default;
        RationalReciprocalLatticeDirection(RationalReciprocalLatticeDirection<dim>&& other) =default;
        VectorDimD cartesian() const;
        Rational<IntScalarType> dot(const LatticeVector<dim>& other) const;
        RationalReciprocalLatticeDirection<dim> operator*(const IntScalarType& scalar) const;
        RationalReciprocalLatticeDirection<dim> operator/(const IntScalarType& scalar) const;
        RationalReciprocalLatticeDirection<dim> operator+(const RationalReciprocalLatticeDirection<dim>& other) const;
        RationalReciprocalLatticeDirection<dim> operator-(const RationalReciprocalLatticeDirection<dim>& other) const;
        RationalReciprocalLatticeDirection<dim> operator+(const ReciprocalLatticeVector<dim>& other) const;
        RationalReciprocalLatticeDirection<dim> operator-(const ReciprocalLatticeVector<dim>& other) const;
        double squaredNorm() const;
        
    };
    
    template<int dim>
    RationalReciprocalLatticeDirection<dim> operator*(const typename RationalReciprocalLatticeDirection<dim>::IntScalarType& scalar, const RationalReciprocalLatticeDirection<dim>& L);
    
} // end namespace
#endif
