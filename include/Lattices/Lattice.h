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
#include <range.h>
#include <map>
#include <RLLL.h>

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
        std::vector<LatticeDirection<dim>> planeParallelLatticeBasis(const ReciprocalLatticeDirection<dim>& l,
                                                                     const bool& useRLLL=false) const;


        /*! \brief Given a lattice direction \f$\textbf l\f$, this function returns a direction-orthogonal reciprocal
         * lattice basis \f$[\textbf r_1,\cdots,\textbf r_{dim}]\f$, with the property
         *  \f$\textbf r_1 \cdot \textbf l=1\f$ and the remaining reciprocal basis vectors are orthogonal to \f$\textbf l\f$, i.e.,
         *  \f$\boldsymbol r_i \cdot \textbf l = 0\f$ for \f$i=2,\cdots,dim\f$.
         *
         * \param[in] l Lattice direction
         *  \returns  a reciprocal lattice basis \f$[\textbf r_1,\cdots,\textbf r_{dim}]\f$
         * */
        std::vector<ReciprocalLatticeDirection<dim>> directionOrthogonalReciprocalLatticeBasis(const LatticeDirection<dim>& l,
                                                                                               const bool& useRLLL=false) const;

        ReciprocalLatticeDirection<dim> reciprocalLatticeDirection(const VectorDimD& d) const;
        RationalLatticeDirection<dim> rationalLatticeDirection(const VectorDimD& d,
                                                               const typename BestRationalApproximation::LongIntType& maxDen=1000) const;
        RationalReciprocalLatticeDirection<dim> rationalReciprocalLatticeDirection(const VectorDimD& d,
                                                               const typename BestRationalApproximation::LongIntType& maxDen=1000) const;

        LatticeVector<dim> latticeVector(const VectorDimD& p) const;
        ReciprocalLatticeVector<dim> reciprocalLatticeVector(const VectorDimD& p) const;

        /*!
         * \brief Computes the interplaar spacing
         * @param r Reciprocal lattice direction
         * @return interPlanarSpacing
         */
        double interPlanarSpacing(const ReciprocalLatticeDirection<dim>& r) const;

        /*! This function generates lattices that share a coincidence relation with the current
         * lattice using rotations about a given axis. It is specialized to dim=3
         *
         * @tparam dm dimension (int)
         * @param rd axis (Reciprocal lattice direction)
         * @param maxStrain
         * @return A set of rotated lattices that are coincident with the current lattice.
         */
        template<int dm=dim>
        typename std::enable_if<dm==3,std::map<IntScalarType ,Lattice<dm>>>::type
        generateCoincidentLattices(const ReciprocalLatticeDirection<dim>& rd, const double& maxStrain=0.0) const
        {
            std::map<IntScalarType,Lattice<dim>> output;
            auto basis= planeParallelLatticeBasis(rd);
            const int N=100;
            const int maxDen=100;
            double epsilon=1e-8;
            IntScalarType keyScale= 1e6;

            auto b1= basis[1].cartesian();
            auto b2= basis[2].cartesian();
            Eigen::Matrix<double,3,2> b;
            b.col(0)= b1; b.col(1)= b2;
            //RLLL(b,0.75).reducedBasis();
            //b1= b.col(0);
            //b2= b.col(1);

            for(IntScalarType i : range<IntScalarType>(-N,N))
            {
                for(IntScalarType j : range<IntScalarType>(-N,N))
                {
                    auto vec = i*b1+j*b2;
                    if (abs(vec.norm())<epsilon ) continue;
                    double ratio= vec.norm()/b1.norm();
                    BestRationalApproximation bra(ratio,maxDen);
                    double error= ratio-static_cast<double>(bra.num)/bra.den;
                    if (abs(error) > epsilon)
                        continue;
                    else
                    {
                        double cosTheta= b1.dot(vec)/(b1.norm()*vec.norm());
                        if (cosTheta-1>-epsilon) cosTheta= 1.0;
                        if (cosTheta+1<epsilon) cosTheta= -1.0;
                        double theta= acos(cosTheta);
                        Eigen::Matrix3d rotation;
                        rotation= Eigen::AngleAxis<double>(theta,rd.cartesian().normalized());
                        IntScalarType key= theta*keyScale;
                        output.insert(std::pair<IntScalarType,Lattice<dim>>(key,Lattice<dim>(this->latticeBasis,rotation)));
                    }
                }
            }
            return output;
        }

        template<int dm=dim>
        typename std::enable_if<dm==2,std::vector<Lattice<dm>>>::type generateCoincidentLattices(const double& maxStrain=0.0) const
        {
            std::vector<Lattice<dim>> output;
            return output;
        }
};
/*! @example testPlaneParallelLatticeDirections.cpp
 *  This example demonstrates the computation of plane-parallel lattice basis and direction-orthogonal reciprocal
 *  lattice basis for a random reciprocal lattice direction and a random lattice direction, respectively.
 *
 *  -# Set up the random number distribution for generating random input reciprocal and lattice directions
 *  @snippet testPlaneParallelLatticeDirections.cpp Random
 *
 *  -# Instantiate a lattice
 *  @snippet testPlaneParallelLatticeDirections.cpp Lattice
 *
 *  -# Form random reciprocal lattice and lattice directions
 *  @snippet testPlaneParallelLatticeDirections.cpp Directions
 *
 *  -# Compute the plane-parallel lattice basis
 *  @snippet testPlaneParallelLatticeDirections.cpp Basis1
 *
 *  -# Test the plane-parallel lattice basis
 *  @snippet testPlaneParallelLatticeDirections.cpp Test1
 *
 *  -# Compute the direction-orthogonal reciprocal lattice basis
 *  @snippet testPlaneParallelLatticeDirections.cpp Basis2
 *
 *  -# Test the direction-orthogonal reciprocal lattice basis
 *  @snippet testPlaneParallelLatticeDirections.cpp Test2
 *
 *
 *  Full Code:
*/

/*! @example testCoincidentRotations.cpp
 *  This example demonstrates the generation of CSLs using rotations about an axis of a 3D lattice
 *
 * -# Define types
 * @snippet testCoincidentRotations.cpp Types
 *
 * -# Instantiate a lattice \f$\mathcal A\f$
 * @snippet testCoincidentRotations.cpp Lattice
 *
 * -# Specify an axis in the form of a reciprocal lattice direction
 * @snippet testCoincidentRotations.cpp Axis
 *
 * -# Generate all rotations \f$\mathbf R\f$ about the given axis that result in a coincidence relation between
 * \f$\mathcal A\f$ and \f$\mathbf R\mathcal A\f$
 * @snippet testCoincidentRotations.cpp Test
 *
 * -# SNF bicrystallography of the resulting bicrystals
 * @snippet testCoincidentRotations.cpp SNF
 *
 *
 * Full code:
 */

/*!
 * @example testLattice.cpp
 * This example demonstrates the use of Lattice class
 *
 * -# Initializing a lattice
 *  @snippet testLattice.cpp Lattice
 *
 * -# Initializing a lattice vector
 *  @snippet testLattice.cpp Lattice vector
 *
 *  -# Modifying its integer coordinates
 *   @snippet testLattice.cpp Integer coordinates
 *
 *  -# Initializing using Cartesian coordinates
 *   @snippet testLattice.cpp Cartesian coordinates
 *
 * -# Lattice vector algebra
 *  @snippet testLattice.cpp Lattice vector algebra
 *
 *
 * Full code:
 */
}
#endif
