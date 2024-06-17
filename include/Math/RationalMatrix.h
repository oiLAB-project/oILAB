/* This file is part of gbLAB.
 *
 * gbLAB is distributed without any warranty under the MIT License.
 */


#ifndef gbLAB_RationalMatrix_h_
#define gbLAB_RationalMatrix_h_

#include <Eigen/Dense>
namespace gbLAB
{
    
    template <int dim>
    class RationalMatrix
    {

        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        typedef long long int IntScalarType;
        typedef Eigen::Matrix<IntScalarType,dim,1> VectorDimI;
        typedef Eigen::Matrix<IntScalarType,dim,dim> MatrixDimI;

        //static constexpr int64_t maxDen=10000000;
        static constexpr long long int maxDen=1000000;
        static std::pair<MatrixDimI,IntScalarType> compute(const MatrixDimD& R);
        static std::pair<MatrixDimI,IntScalarType> reduce(const MatrixDimI& Rn, const MatrixDimI& Rd);

        const std::pair<MatrixDimI,IntScalarType> returnPair;
        
    public:
        
        const MatrixDimI& integerMatrix;
        const IntScalarType& mu;
        
        RationalMatrix(const MatrixDimD& R) ;
        RationalMatrix(const MatrixDimI& Rn,const IntScalarType& Rd);
        RationalMatrix(const MatrixDimI& Rn, const MatrixDimI& Rd);
        MatrixDimD asMatrix() const;

    };

} // end namespace
#endif


