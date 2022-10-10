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

        /*! \brief Obtain plane-parallel lattice directions
         *
         * \param[in] i Miller indices of a lattice plane
         *  \returns  dim-1 planar lattice directions
         *  \f$x=\sqrt{y}\f$
         * */
        std::vector<LatticeDirection<dim>> planeParallelLatticeBasis(const ReciprocalLatticeDirection<dim>& l) const;

        /*! \brief Obtain plane-parallel lattice directions
         *
         * \param[in] i integer coordinates of a lattice direction
         *  \returns  dim-1 reciprocal lattice directions perpendicular to the input lattice direction
         *  \f$x=\sqrt{y}\f$
         * */
        std::vector<ReciprocalLatticeDirection<dim>> directionOrthogonalReciprocalLatticeBasis(const LatticeDirection<dim>& l) const;

        ReciprocalLatticeDirection<dim> reciprocalLatticeDirection(const VectorDimD& d) const;
        RationalLatticeDirection<dim> rationalLatticeDirection(const VectorDimD& d,
                                                              const typename BestRationalApproximation::LongIntType& maxDen=1000) const;

        LatticeVector<dim> latticeVector(const VectorDimD& p) const;
        ReciprocalLatticeVector<dim> reciprocalLatticeVector(const VectorDimD& p) const;

};
}
#endif
