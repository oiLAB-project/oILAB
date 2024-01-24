/* This file is part of gbLAB.
 *
 * gbLAB is distributed without any warranty under the MIT License.
 */

#ifndef gbLAB_LatticeVector_h_
#define gbLAB_LatticeVector_h_

#include <LatticeModule.h>

namespace gbLAB
{
    /*! \brief LatticeVector class
     *
     *  The LatticeVector<dim> class describes a lattice vector in a lattice
     * */
    template <int dim>
    class LatticeVector: public Eigen::Matrix<typename LatticeCore<dim>::IntScalarType,dim,1>
    {
        typedef Eigen::Matrix<typename LatticeCore<dim>::IntScalarType,dim,1> BaseType;
        BaseType& base();
        const BaseType& base() const;

    public:
        
//        static constexpr double roundTol=FLT_EPSILON;
        typedef typename LatticeCore<dim>::IntScalarType IntScalarType;
        typedef typename LatticeCore<dim>::VectorDimD VectorDimD;
        typedef typename LatticeCore<dim>::MatrixDimD MatrixDimD;
        typedef typename LatticeCore<dim>::VectorDimI VectorDimI;
        typedef typename LatticeCore<dim>::MatrixDimI MatrixDimI;
        
        const Lattice<dim>& lattice;
        
        LatticeVector(const Lattice<dim>& lat) ;
        LatticeVector(const VectorDimD& d, const Lattice<dim>& lat) ;
        LatticeVector(const VectorDimI& d, const Lattice<dim>& lat) ;
        /*
        template<typename OtherDerived>
        LatticeVector(const Eigen::MatrixBase<OtherDerived>& other, const Lattice<dim>& lat) :
        BaseType(other),
        lattice(lat)
        {
        }
        */
        LatticeVector(const LatticeVector<dim>& other) = default;
        LatticeVector(LatticeVector<dim>&& other) =default;
        LatticeVector<dim>& operator=(const LatticeVector<dim>& other);
        LatticeVector<dim>& operator=(LatticeVector<dim>&& other);
        LatticeVector<dim> operator+(const LatticeVector<dim>& other) const;
        LatticeVector<dim>& operator+=(const LatticeVector<dim>& other);
        LatticeVector<dim> operator-(const LatticeVector<dim>& other) const;
        LatticeVector<dim>& operator-=(const LatticeVector<dim>& other);
        LatticeVector<dim> operator*(const IntScalarType& scalar) const
        {
            VectorDimI temp= static_cast<VectorDimI>(*this) * scalar;
            return LatticeVector<dim>(temp, lattice);
        }
        IntScalarType dot(const ReciprocalLatticeVector<dim>& other) const;
        IntScalarType dot(const ReciprocalLatticeDirection<dim>& other) const;
        VectorDimD cartesian() const;

        template<int dm=dim>
        typename std::enable_if<dm==3,ReciprocalLatticeDirection<dm>>::type
        cross(const LatticeVector<dm>& other) const
        {
            assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
            return ReciprocalLatticeDirection<dm>(ReciprocalLatticeVector<dm>(static_cast<VectorDimI>(*this).cross(static_cast<VectorDimI>(other)), lattice));
        }

    };
        
    template<int dim>
    LatticeVector<dim> operator*(const typename LatticeVector<dim>::IntScalarType& scalar, const LatticeVector<dim>& L);

    template<int dim>
    LatticeVector<dim> operator*(const int& scalar, const LatticeVector<dim>& L);
} // end namespace
#endif
