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
    template<int dim>
    LatticeVector<dim> LatticeVector<dim>::operator*(const LatticeVector<dim>::IntScalarType& scalar) const
    {
        VectorDimI temp= static_cast<VectorDimI>(*this) * scalar;
        return LatticeVector<dim>(temp, lattice);
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
    typename LatticeVector<dim>::IntScalarType LatticeVector<dim>::dot(const ReciprocalLatticeDirection<dim> &other) const
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        return dot(other.reciprocalLatticeVector());
    }

    /**********************************************************************/
    template <int dim>
    typename LatticeVector<dim>::VectorDimD LatticeVector<dim>::cartesian() const
    {
        return lattice.latticeBasis * this->template cast<double>();
    }

    /**********************************************************************/
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


    template<int dim> template<int dm>
    typename std::enable_if<dm==3,void>::type
    LatticeVector<dim>::modulo(LatticeVector<dim>& input, const std::vector<LatticeVector<dim>>& basis, const VectorDimD& shift)
    {
        double det= (basis[0].cross(basis[1])).dot(basis[2]);
        assert( abs(det) > FLT_EPSILON  );
        auto normal= basis[1].cross(basis[2]);
        input= input - floor( (double)input.dot(normal)/basis[0].dot(normal)-shift(0) ) * basis[0];
        assert((double)(input.dot(normal))/basis[0].dot(normal)<= shift(0)+1.0 &&
               (double)(input.dot(normal))/basis[0].dot(normal)>= shift(0));
        normal= basis[2].cross(basis[0]);
        input= input - floor( (double)input.dot(normal)/basis[1].dot(normal)-shift(1) ) * basis[1];
        assert((double)(input.dot(normal))/basis[1].dot(normal)<= shift(1)+1.0 &&
               (double)(input.dot(normal))/basis[1].dot(normal)>= shift(1));
        normal= basis[0].cross(basis[1]);
        input= input - floor( (double)input.dot(normal)/basis[2].dot(normal)-shift(2) ) * basis[2];
        assert((double)(input.dot(normal))/basis[2].dot(normal)<= shift(2)+1.0 &&
               (double)(input.dot(normal))/basis[2].dot(normal)>= shift(2));
    }

    template<int dim> template<int dm>
    typename std::enable_if<dm==3,void>::type
    LatticeVector<dim>::modulo(VectorDimD& input, const std::vector<LatticeVector<dim>>& basis, const VectorDimD& shift)
    {
        Eigen::Matrix3d L;
        L.col(0)= basis[0].cartesian();
        L.col(1)= basis[1].cartesian();
        L.col(2)= basis[2].cartesian();

        Eigen::Vector3d inputCoordinates= ((L.inverse()*input).array()-shift.array()).floor();
        input= input - L*inputCoordinates;
    }


    template<int dim> template<int dm>
    typename std::enable_if<dm==2,void>::type
    LatticeVector<dim>::modulo(LatticeVector<dim>& input, const std::vector<LatticeVector<dim>>& basis, const VectorDimD& shift)
    {
        auto normal= basis[1].cross();
        input= input - input.dot(normal)/basis[0].dot(normal) * basis[0];
        normal= basis[0].cross();
        input= input - input.dot(normal)/basis[1].dot(normal) * basis[1];
    }

    template<int dim> template<int dm>
    typename std::enable_if<dm==2,void>::type
    LatticeVector<dim>::modulo(VectorDimD& input, const std::vector<LatticeVector<dim>>& basis, const VectorDimD& shift)
    {
    }


    template class LatticeVector<1>;
    template LatticeVector<1>operator*(const typename LatticeVector<1>::IntScalarType& scalar, const LatticeVector<1>& L);
    template LatticeVector<1>operator*(const int& scalar, const LatticeVector<1>& L);

    template class LatticeVector<2>;
    template LatticeVector<2>operator*(const typename LatticeVector<2>::IntScalarType& scalar, const LatticeVector<2>& L);
    template LatticeVector<2>operator*(const int& scalar, const LatticeVector<2>& L);
    template void LatticeVector<2>::modulo<2>(LatticeVector<2>& input, const std::vector<LatticeVector<2>>& basis, const Eigen::Vector2d& shift);
    template void LatticeVector<2>::modulo<2>(Eigen::Vector2d& input, const std::vector<LatticeVector<2>>& basis, const Eigen::Vector2d& shift);

    template class LatticeVector<3>;
    template LatticeVector<3>operator*(const typename LatticeVector<3>::IntScalarType& scalar, const LatticeVector<3>& L);
    template LatticeVector<3>operator*(const int& scalar, const LatticeVector<3>& L);
    template void LatticeVector<3>::modulo<3>(LatticeVector<3>& input, const std::vector<LatticeVector<3>>& basis, const Eigen::Vector3d& shift);
    template void LatticeVector<3>::modulo<3>(Eigen::Vector3d& input, const std::vector<LatticeVector<3>>& basis, const Eigen::Vector3d& shift);

    template class LatticeVector<4>;
    template LatticeVector<4>operator*(const typename LatticeVector<4>::IntScalarType& scalar, const LatticeVector<4>& L);
    template LatticeVector<4>operator*(const int& scalar, const LatticeVector<4>& L);

    template class LatticeVector<5>;
    template LatticeVector<5>operator*(const typename LatticeVector<5>::IntScalarType& scalar, const LatticeVector<5>& L);
    template LatticeVector<5>operator*(const int& scalar, const LatticeVector<5>& L);
} // end namespace
#endif
