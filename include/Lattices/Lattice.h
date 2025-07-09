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
#include <unordered_map>
#include "RLLL.h"
#include <Rational.h>
#include <Farey.h>
#include <RationalApproximations.h>
#include <algorithm>
#include <fstream>


namespace gbLAB
{
    /*! \brief Lattice class
     *
     *  The Lattice<dim> class describes a lattice in dim dimensions
     * */
    template <int dim>
    class Lattice : public StaticID<Lattice<dim>>
    {
        static constexpr double roundTol=FLT_EPSILON;

        using VectorDimD = typename LatticeCore<dim>::VectorDimD;
        using MatrixDimD = typename LatticeCore<dim>::MatrixDimD;
        using VectorDimI = typename LatticeCore<dim>::VectorDimI;
        using MatrixDimI = typename LatticeCore<dim>::MatrixDimI;
        using IntScalarType = typename LatticeCore<dim>::IntScalarType;


    public:
        
        const MatrixDimD    latticeBasis;
        const MatrixDimD reciprocalBasis;
        const MatrixDimD F;

        Lattice(const MatrixDimD& A,const MatrixDimD& Q=MatrixDimD::Identity()) ;


        /*! \brief Returns a lattice vector (in the current lattice) with Cartesian coordinates p
         *
         * @param[in] d cartesian coordinates of a vector
         * @return Lattice vector
         */
        LatticeVector<dim> latticeVector(const VectorDimD& p) const;

        /*! \brief Returns the lattice direction along a vector
         *
         * @param[in] d cartesian coordinates of a vector
         * @return Lattice direction along d
         */
        LatticeDirection<dim> latticeDirection(const VectorDimD& d) const;

        /*! \brief Returns a reciprocal lattice vector (in the dual of the current lattice) with Cartesian coordinates p
         *
         * @param[in] p cartesian coordinates of a vector
         * @return Reciprocal lattice vector
         */
        ReciprocalLatticeVector<dim> reciprocalLatticeVector(const VectorDimD& p) const;

        /*! \brief Returns the reciprocal lattice direction along a vector
         *
         * @param[in] d cartesian coordinates of a vector
         * @return Reciprocal lattice direction along d
         */
        ReciprocalLatticeDirection<dim> reciprocalLatticeDirection(const VectorDimD& d) const;

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

        RationalLatticeDirection<dim> rationalLatticeDirection(const VectorDimD& d,
                                                               const typename BestRationalApproximation::LongIntType& maxDen=1000) const;
        RationalReciprocalLatticeDirection<dim> rationalReciprocalLatticeDirection(const VectorDimD& d,
                                                                                   const typename BestRationalApproximation::LongIntType& maxDen=1000) const;

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



        /*!
         * \brief Computes the interplanar spacing
         * @param r Reciprocal lattice direction
         * @return interPlanarSpacing
         */
        double interPlanarSpacing(const ReciprocalLatticeDirection<dim>& r) const;

        /*! This function generates rotations (about a given axis) that result in a coincident site lattice.
         *  It is specialized to dim=3
         *
         * @tparam dm dimension (int)
         * @param rd axis (Reciprocal lattice direction)
         * @param maxDen  integer parameter that determines the resolution for the search of rotations
         * @param N integer parameter that determines the maximum size of the CSL
         * @return A set of rotations that result in CSLs.
         */
        template<int dm=dim>
        typename std::enable_if<dm==3,std::vector<MatrixDimD>>::type
        generateCoincidentLattices(const ReciprocalLatticeDirection<dim>& rd, const double& maxDen= 100, const int& N= 100) const;

        /*! This function generates deformations \f$\mathbf F\f$ such that the deformations of *this lattice share moire supercells
         *  with the undeformed *this lattice
         */
        template<int dm=dim>
        typename std::enable_if<dm==2,std::vector<MatrixDimD>>::type
        generateCoincidentLattices(const double& maxStrain,
                                   const int& maxDen=50,
                                   const int& N=30) const;

        /*! This function generates deformations \f$\mathbf F\f$ such that the deformations of *this lattice share moire supercells
         *  with a given undeformed 2D lattice. It is specialized to dim=2.
         *
         *  The current algorithm betters the one given in Algorithm 2 of
         *
         *  [1] Admal, Nikhil Chandra, et al. "Interface dislocations and grain boundary
         *      disconnections using Smith normal bicrystallography."
         *      Acta materialia 240 (2022): 118340.
         *
         *  Note: There was a typo in Algorithm 2 in [1] - \f$\mathfrak q_1\f$ and \f$\mathfrak r_1\f$ should
         *  be replaced by the basis vectors \f$\mathbf q_1\f$ and \f$\mathbf r_1\f$ of lattice \f$\mathcal B\f$.
         *
         * @tparam undeformedLattice underformed lattice
         * @tparam dm dimension (int)
         * @param maxStrain maximum strain
         * @param maxDen  integer parameter that determines the resolution for the search of rotations
         * @param N integer parameter that determines the maximum size of the CSL
         * @return A set of deformation gradients of *this lattice that result in moire superlattices with the undeformed
         *         lattice
         */
        template<int dm=dim>
        typename std::enable_if<dm==2,std::vector<MatrixDimD>>::type
        generateCoincidentLattices(const Lattice<dim>& undeformedLattice,
                                   const double& maxStrain,
                                   const int& maxDen=50,
                                   const int& N=30) const;

        /*! This function outputs/prints lattice points within a box bounded by the
         * input box vectors. The box vectors have to be linearly independent lattice
         * vectors. This function is specialized to dim=3.
         *
         * @tparam dm dimension (int)
         * @param boxVectors three linearly independent lattice vectors
         * @param filename (optional) name of the output file
         * @return Lattice points bounded by the box vectors
         */
        template<int dm=dim>
        typename std::enable_if<dm==3,std::vector<LatticeVector<dim>>>::type
        box(const std::vector<LatticeVector<dim>>& boxVectors, const std::string& filename= "") const;

        /*! This function outputs/prints lattice points within a box bounded by the
         * optional input box vectors. The box vectors have to be linearly independent lattice
         * vectors. This function is specialized to dim=2.
         *
         * @tparam dm dimension (int)
         * @param boxVectors two linearly independent lattice vectors
         * @param filename (optional) name of the output file
         * @return Lattice points bounded by the box vectors
         */
        template<int dm=dim>
        typename std::enable_if<dm==2,std::vector<LatticeVector<dim>>>::type
        box(const std::vector<LatticeVector<dim>>& boxVectors, const std::string& filename= "") const;
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
 * -# Generate a family of rotations about the given axis such that lattices \f$\mathbf R\mathcal A\f$ and
 * \f$\mathbf R\mathcal A\f$ share a coincidence relation
 * @snippet testCoincidentRotations.cpp Test
 *
 * -# Loop over the the above family of rotations and carry out SNF bicrystallography of bicrystals
 * \f$\mathcal A \cup \mathbf R\mathcal A\f$.
 * @snippet testCoincidentRotations.cpp SNF
 *
 *
 * Full code:
 */

/*! @example testMoire.cpp
 *  This example demonstrates the generation of strained moire superlattices from a 2D homostructure and calculation
 *  of the translational invariance of a moire
 *
 * -# Define types
 * @snippet testMoire.cpp Types
 *
 * -# Instantiate a lattice \f$\mathcal A\f$
 * @snippet testMoire.cpp Lattice
 *
 * -# Generate deformation gradients \f$\mathbf F\f$ (with Lagrangian strain < 0.01) , such that lattices \f$\mathcal A\f$
 * and \f$\mathbf F\mathcal A\f$ form a moire superlattice
 * @snippet testMoire.cpp Test
 *
 * -# Loop over the deformations to form 2D homostructures, \f$\mathcal A \cup \mathbf F\mathcal A\f$, and carry
 * out SNF bicrystallography
 * @snippet testMoire.cpp SNF
 *
 * -# Output the heterodeformation, its polar decomposition, and the corresponding elastic strain
 * @snippet testMoire.cpp Heterodeform
 *
 * -# Output the invariance property of the moire in the following steps:
 * 1) compute reduced basis vectors \f$\mathbf d_1\f$ and \f$\mathbf d_2\f$ for the DSCL,
 * 2) compute the moire shifts \f$\mathbf s_1\f$ and \f$\mathbf s_2\f$ when lattice \f$\mathcal A\f$ is
 * displaced by \f$\mathbf d_1\f$ and \f$\mathbf d_2\f$, respectively.
 * @snippet testMoire.cpp Invariance
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
 *  -# Initializing using Cartesian coordinates may fail if they don't describe a lattice point
 *   @snippet testLattice.cpp Cartesian coordinates fail
 *
 * -# Lattice vector algebra
 *  @snippet testLattice.cpp Lattice vector algebra
 *
 * -# Lattice direction from a lattice vector
 *  @snippet testLattice.cpp Lattice direction
 *
 * -# Creating a reciprocal lattice vector
 *  @snippet testLattice.cpp Reciprocal lattice vector
 *
 * -# Reciprocal lattice direction from a reciprocal lattice vector
 *  @snippet testLattice.cpp Reciprocal lattice direction
 *
 * -# Get reciprocal lattice vector from a reciprocal lattice direction. This
 * may be required to access functions that accept a vector as an input.
 *  @snippet testLattice.cpp Direction to vector
 *
 * -# Get a lattice direction along a given cartesian vector.
 *  CAUTION: This will fail if the cartesian vector is not a lattice vector.
 *  @snippet testLattice.cpp Lattice direction along a cartesian vector
 *
 * -# Get the stacking along a lattice plane
 *  @snippet testLattice.cpp Stacking of a lattice plane

 * -# The cross product of two lattice vectors is a reciprocal lattice direction.
 * Similarly, the cross product of two reciprocal lattice vectors is a lattice direction.
 * The cross product is enabled only for dim=3
 *  @snippet testLattice.cpp Cross product
 *
 * -# Output a configuration of lattice points bounded by three lattice vectors (named boxVectors)
 * of a 3D lattice.
 *  @snippet testLattice.cpp Box3
 *
 * -# Create a 2D lattice and output its lattice points bounded by two lattice vectors (named boxVectors)
 * of a 2D lattice.
 *  @snippet testLattice.cpp Box2
 *
 * Full code:
 */
}
#endif
