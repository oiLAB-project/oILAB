/* This file is part of gbLAB.
 *
 * gbLAB is distributed without any warranty under the MIT License.
 */


#ifndef gbLAB_Lattice_h_
#define gbLAB_Lattice_h_

#include <LatticeModule.h>
#include <StaticID.h>
#include <BestRationalApproximation.h>
#include <vector>
#include <IntegerLattice.h>

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
        template<int d=dim>
        typename std::enable_if< (d>1 && d<4),std::vector<LatticeDirection<d>>> :: type
        planeParallelLatticeDirections(const ReciprocalLatticeDirection<d>& l) const
        {
            auto intLatticeDirections= IntegerLattice<d>::perpendicularDirections(l);
            std::vector<LatticeDirection<d>> out;
            for(int ind= 0; ind<d-1; ind++)
            {
                out.push_back(LatticeDirection<d>(intLatticeDirections.row(ind),*this));
            }
            return out;
        }

        /*! \brief Obtain plane-parallel lattice directions
         *
         * \param[in] i integer coordinates of a lattice direction
         *  \returns  dim-1 reciprocal lattice directions perpendicular to the input lattice direction
         *  \f$x=\sqrt{y}\f$
         * */
        template<int d=dim>
        typename std::enable_if< (d>1 && d<4),std::vector<ReciprocalLatticeDirection<d>>> :: type
        planeParallelLatticeDirections(const LatticeDirection<d>& l) const
        {
            auto intReciprocalLatticeDirections= IntegerLattice<d>::perpendicularDirections(l);
            std::vector<ReciprocalLatticeDirection<d>> out;
            for(int ind= 0; ind<d-1; ind++)
            {
                out.push_back(ReciprocalLatticeDirection<d>(intReciprocalLatticeDirections.row(ind),*this));
            }
            return out;
        }

        template<typename OtherDerived>
        LatticeDirection<dim> latticeDirection(const Eigen::MatrixBase<OtherDerived>& other) const
        {
            return LatticeDirection<dim>(LatticeVector<dim>(other,*this));
        }
        
        ReciprocalLatticeDirection<dim> reciprocalLatticeDirection(const VectorDimD& d) const;
        RationalLatticeDirection<dim> rationalLatticeDirection(const VectorDimD& d,
                                                              const typename BestRationalApproximation::LongIntType& maxDen=1000) const;
        
        template<typename OtherDerived>
        LatticeVector<dim> latticeVector(const Eigen::MatrixBase<OtherDerived>& other) const
        {
            return LatticeVector<dim>(other,*this);
        }
        
        LatticeVector<dim> latticeVector(const VectorDimD& p) const;
        ReciprocalLatticeVector<dim> reciprocalLatticeVector(const VectorDimD& p) const;

    };
}
#endif
