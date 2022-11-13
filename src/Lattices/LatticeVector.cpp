/* This file is part of gbLAB.
 *
 * gbLAB is distributed without any warranty under the MIT License.
 */


#ifndef gbLAB_LatticeVector_cpp_
#define gbLAB_LatticeVector_cpp_

#include<LatticeModule.h>

namespace gbLAB
{

    /**********************************************************************/
    template <int dim>
    typename LatticeVector<dim>::BaseType& LatticeVector<dim>::base()
    {
        return *this;
    }

    /**********************************************************************/
    template <int dim>
    const typename LatticeVector<dim>::BaseType& LatticeVector<dim>::base() const
    {
        return *this;
    }


    /**********************************************************************/
    template <int dim>
    LatticeVector<dim>::LatticeVector(const Lattice<dim> &lat) :
    /* init */ BaseType(VectorDimI::Zero()),
    /* init */ lattice(lat)
    {
    }

    /**********************************************************************/
    template <int dim>
    LatticeVector<dim>::LatticeVector(const VectorDimD &d,
                  const Lattice<dim> &lat) :
    /* init */ BaseType(LatticeCore<dim>::integerCoordinates(d,lat.reciprocalBasis.transpose())),
    /* init */ lattice(lat)
    { /*!@param[in] d vector in real space
          * Constructs *this by mapping d to the lattice
          */
    }

    /**********************************************************************/
    template<int dim>
    LatticeVector<dim>::LatticeVector(const VectorDimI& other, const Lattice<dim>& lat) :
            /* init base */ BaseType(other),
            /* init      */ lattice(lat)
    { }

    /**********************************************************************/
    template <int dim>
    LatticeVector<dim>& LatticeVector<dim>::operator=(const LatticeVector<dim> &other)
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        base() = other.base();
        return *this;
    }
    /**********************************************************************/
    template <int dim>
    LatticeVector<dim>& LatticeVector<dim>::operator=(LatticeVector<dim> &&other)
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        base() = other.base();
        return *this;
    }

    /**********************************************************************/
    template <int dim>
    LatticeVector<dim> LatticeVector<dim>::operator+(const LatticeVector<dim> &other) const
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        VectorDimI temp= static_cast<VectorDimI>(*this) + static_cast<VectorDimI>(other);
        return LatticeVector<dim>(temp, lattice);
    }

    /**********************************************************************/
    template <int dim>
    LatticeVector<dim>& LatticeVector<dim>::operator+=(const LatticeVector<dim> &other)
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        base() += other.base();
        return *this;
    }

    /**********************************************************************/
    template <int dim>
    LatticeVector<dim> LatticeVector<dim>::operator-(const LatticeVector<dim> &other) const
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        VectorDimI temp= static_cast<VectorDimI>(*this) - static_cast<VectorDimI>(other);
        return LatticeVector<dim>(temp, lattice);
    }

    /**********************************************************************/
    template <int dim>
    LatticeVector<dim>& LatticeVector<dim>::operator-=(const LatticeVector<dim> &other)
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        base() -= other.base();
        return *this;
    }

    /**********************************************************************/
    template <int dim>
    typename LatticeVector<dim>::IntScalarType LatticeVector<dim>::dot(const ReciprocalLatticeVector<dim> &other) const
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        return static_cast<VectorDimI>(*this).dot(static_cast<VectorDimI>(other));
    }

    /**********************************************************************/
    template <int dim>
    typename LatticeVector<dim>::VectorDimD LatticeVector<dim>::cartesian() const
    {
        return lattice.latticeBasis * this->template cast<double>();
    }
    
    template<int dim>
    LatticeVector<dim> operator*(const typename LatticeVector<dim>::IntScalarType& scalar, const LatticeVector<dim>& L)
    {
        return L*scalar;
    }

    template<int dim>
    LatticeVector<dim> operator*(const int& scalar, const LatticeVector<dim>& L)
    {
        return L*scalar;
    }

    template class LatticeVector<1>;
    template LatticeVector<1>operator*(const typename LatticeVector<1>::IntScalarType& scalar, const LatticeVector<1>& L);
    template LatticeVector<1>operator*(const int& scalar, const LatticeVector<1>& L);
    template class LatticeVector<2>;
    template LatticeVector<2>operator*(const typename LatticeVector<2>::IntScalarType& scalar, const LatticeVector<2>& L);
    template LatticeVector<2>operator*(const int& scalar, const LatticeVector<2>& L);
    template class LatticeVector<3>;
    template LatticeVector<3>operator*(const typename LatticeVector<3>::IntScalarType& scalar, const LatticeVector<3>& L);
    template LatticeVector<3>operator*(const int& scalar, const LatticeVector<3>& L);
    template class LatticeVector<4>;
    template LatticeVector<4>operator*(const typename LatticeVector<4>::IntScalarType& scalar, const LatticeVector<4>& L);
    template LatticeVector<4>operator*(const int& scalar, const LatticeVector<4>& L);
    template class LatticeVector<5>;
    template LatticeVector<5>operator*(const typename LatticeVector<5>::IntScalarType& scalar, const LatticeVector<5>& L);
    template LatticeVector<5>operator*(const int& scalar, const LatticeVector<5>& L);
} // end namespace
#endif
