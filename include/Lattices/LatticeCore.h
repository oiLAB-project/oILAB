/* This file is part of gbLAB.
 *
 * gbLAB is distributed without any warranty under the MIT License.
 */


#ifndef gbLAB_LatticeCore_h_
#define gbLAB_LatticeCore_h_

#include <iostream>
#include <cfloat> // FLT_EPSILON
#include <Eigen/Dense>
#include <IntegerMath.h>


namespace gbLAB
{
    template <int dim>
    struct LatticeCore
    {
        static_assert(dim>0,"dim must be > 0.");
        static constexpr double roundTol=FLT_EPSILON;

        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        typedef long long int IntScalarType;
        typedef Eigen::Matrix<IntScalarType,dim,1> VectorDimI;
        typedef Eigen::Matrix<IntScalarType,dim,dim> MatrixDimI;

        /*!
         * \brief Approximates a direction in terms of integer coordinates
         * @param v input direction
         * @return direction with integer coordinates
         */
        static VectorDimI rationalApproximation(VectorDimD v);

        /*!
         * \brief Returns the integer coordinates of a vector \f$d\f$ with respect to a lattices with
         * structure matrix \f$\textbf A\f$.
         * @param v input direction
         * @param invA \f$\textbf A^{-1}\f$
         * @return integer coordinates
         */
        static VectorDimI integerCoordinates(const VectorDimD& d,const MatrixDimD& invA);
    };
}
#endif
