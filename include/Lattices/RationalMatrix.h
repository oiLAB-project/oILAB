/* This file is part of gbLAB.
 *
 * gbLAB is distributed without any warranty under the MIT License.
 */


#ifndef gbLAB_RationalMatrix_h_
#define gbLAB_RationalMatrix_h_

#include <LatticeCore.h>
#include <BestRationalApproximation.h>

namespace gbLAB
{
    
    template <int dim>
    class RationalMatrix
    {
        typedef typename LatticeCore<dim>::IntScalarType IntScalarType;
        typedef typename LatticeCore<dim>::VectorDimD VectorDimD;
        typedef typename LatticeCore<dim>::MatrixDimD MatrixDimD;
        typedef typename LatticeCore<dim>::VectorDimI VectorDimI;
        typedef typename LatticeCore<dim>::MatrixDimI MatrixDimI;

        static constexpr int64_t maxDen=10000000;
        static std::pair<MatrixDimI,IntScalarType> compute(const MatrixDimD& R);
        
        const std::pair<MatrixDimI,IntScalarType> returnPair;
        
    public:
        
        const MatrixDimI& integerMatrix;
        const IntScalarType& mu;
        
        RationalMatrix(const MatrixDimD& R) ;
        RationalMatrix(const MatrixDimI& Rn,const IntScalarType& Rd);
        MatrixDimD asMatrix() const;

    };

} // end namespace
#endif


