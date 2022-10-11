/* This file is part of gbLAB.
 *
 * gbLAB is distributed without any warranty under the MIT License.
 */


#ifndef gbLAB_Lattice_h_
#define gbLAB_Lattice_h_

#include <IntegerLattice.h>
#include <LatticeModule.h>
#include <StaticID.h>
#include <BestRationalApproximation.h>
#include <vector>

namespace gbLAB
{
    /*! \brief Lattice class
     *
     *  lattice class description
     * */
    template <int dim>
    class Lattice : public StaticID<Lattice<dim>>
    {
        static constexpr double roundTol=FLT_EPSILON;
        typedef typename LatticeCore<dim>::VectorDimD VectorDimD;
        typedef typename LatticeCore<dim>::MatrixDimD MatrixDimD;
        typedef typename LatticeCore<dim>::VectorDimI VectorDimI;
        typedef typename LatticeCore<dim>::MatrixDimI MatrixDimI;
        typedef typename LatticeCore<dim>::IntScalarType IntScalarType;

//        static MatrixDimD getLatticeBasis(const MatrixDimD& A,const MatrixDimD& Q);
        
    public:
        
        const MatrixDimD    latticeBasis;
        const MatrixDimD reciprocalBasis;
        const MatrixDimD F;

        Lattice(const MatrixDimD& A,const MatrixDimD& Q=MatrixDimD::Identity()) ;
        LatticeDirection<dim> latticeDirection(const VectorDimD& d) const;

        /*! \brief Given a reciprocal lattice direction \f$\textbf l\f$, this function returns a plane-parallel lattice basis
         * \f$[\textbf b_1,\cdots,\textbf b_{dim}]\f$, with the property
         *  \f$\textbf b_1 \cdot \textbf l=1\f$ and the remaining basis vectors lie in the lattice plane represented by
         *  \f$\textbf l\f$, i.e. \f$\textbf b_i \cdot \textbf l = 0\f$ for \f$i=2,\cdots,dim\f$.
         *
         * \param[in] l Reciprocal lattice direction
         * \returns   A lattice basis \f$[\textbf b_1,\cdots,\textbf b_{dim}]\f$
         * */
        std::vector<LatticeDirection<dim>> planeParallelLatticeBasis(const ReciprocalLatticeDirection<dim>& l) const;

        /*! \brief Given a lattice direction \f$\textbf l\f$, this function returns a direction-orthogonal reciprocal
         * lattice basis \f$[\textbf r_1,\cdots,\textbf r_{dim}]\f$, with the property
         *  \f$\textbf r_1 \cdot \textbf l=1\f$ and the remaining reciprocal basis vectors are orthogonal to \f$\textbf l\f$, i.e.,
         *  \f$\boldsymbol r_i \cdot \textbf l = 0\f$ for \f$i=2,\cdots,dim\f$.
         *
         * \param[in] l Lattice direction
         *  \returns  a reciprocal lattice basis \f$[\textbf r_1,\cdots,\textbf r_{dim}]\f$
         * */
        std::vector<ReciprocalLatticeDirection<dim>> directionOrthogonalReciprocalLatticeBasis(const LatticeDirection<dim>& l) const;

        ReciprocalLatticeDirection<dim> reciprocalLatticeDirection(const VectorDimD& d) const;
        RationalLatticeDirection<dim> rationalLatticeDirection(const VectorDimD& d,
                                                              const typename BestRationalApproximation::LongIntType& maxDen=1000) const;

        LatticeVector<dim> latticeVector(const VectorDimD& p) const;
        ReciprocalLatticeVector<dim> reciprocalLatticeVector(const VectorDimD& p) const;
//        std::vector<Eigen::Matrix<double,3,3>> transformations(const LatticeDirection<3>& d, const double& maxStrain=0.0) const;

};
}
#endif
