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
#include <RLLL.h>
#include <Rational.h>
#include <Farey.h>
#include <RationalApproximations.h>
#include <algorithm>


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
        generateCoincidentLattices(const ReciprocalLatticeDirection<dim>& rd, const double& maxDen= 100, const int N= 100) const
        {
            std::vector<MatrixDimD> output;
            std::map<IntScalarType,MatrixDimD> temp;
            auto basis= planeParallelLatticeBasis(rd);
            double epsilon=1e-8;
            IntScalarType keyScale= 1e6;

            auto b1= basis[1].cartesian();
            auto b2= basis[2].cartesian();

            // Following the notation in Algorithm 3 of
            // "Interface dislocations and grain boundary disconnections using Smith normal bicrystallography"
            auto q1= b1;
            std::vector<std::pair<int,int>> coPrimePairs= farey(N,false);
            for(const auto& pair : coPrimePairs)
            {
                VectorDimI pIndices;
                auto q2= pair.first*b1+pair.second*b2;

                if (abs(q2.norm())<epsilon ) continue;
                double ratio= q2.norm()/q1.norm();
                BestRationalApproximation bra(ratio,maxDen);
                double error= ratio-static_cast<double>(bra.num)/bra.den;
                if (abs(error) > epsilon)
                    continue;
                else
                {
                    double cosTheta= b1.dot(q2)/(q1.norm()*q2.norm());
                    if (cosTheta-1>-epsilon) cosTheta= 1.0;
                    if (cosTheta+1<epsilon) cosTheta= -1.0;
                    double theta= acos(cosTheta);
                    MatrixDimD rotation;
                    rotation= Eigen::AngleAxis<double>(theta,rd.cartesian().normalized());
                    IntScalarType key= theta*keyScale;
                    temp.insert(std::pair<IntScalarType,MatrixDimD>(key,rotation));
                }
            }
            std::transform(temp.begin(), temp.end(),
                           std::back_inserter(output),
                           [](const std::pair<IntScalarType,MatrixDimD>& p) {
                               return p.second;
                           });
            return output;
        }

        /*! This function generates deformations of a 2D lattice that result in a moire superlattice.
         *  It is specialized to dim=2
         *
         * @tparam dm dimension (int)
         * @param maxStrain maximum strain
         * @param maxDen  integer parameter that determines the resolution for the search of rotations
         * @param N integer parameter that determines the maximum size of the CSL
         * @return A set of deformation gradients that result in moire supercells.
         */
        template<int dm=dim>
        typename std::enable_if<dm==2,std::vector<MatrixDimD>>::type
        generateCoincidentLattices(const double& maxStrain,const int& maxDen=50, const int& N=30) const
        {
            std::vector<MatrixDimD> output;
            std::map<IntScalarType,MatrixDimD> temp;
            MatrixDimI mn(MatrixDimI::Zero());
            MatrixDimI md(MatrixDimI::Ones());;
            std::vector<std::pair<int,int>> coPrimePairs= farey(N,false);
            double epsilon=1e-8;
            IntScalarType keyScale= 1e6;

            for(const auto& pair1 : coPrimePairs)
            {
                VectorDimI pIndices;
                pIndices << pair1.first, pair1.second;
                LatticeVector<dm> q2(pIndices,*this);
                double ratio= q2.cartesian().norm() / latticeBasis.col(0).norm();


                if (abs(maxStrain) < epsilon)
                {
                    BestRationalApproximation alpha(ratio, maxDen);
                    double error = ratio - static_cast<double>(alpha.num) / alpha.den;
                    if (abs(error) > epsilon)
                        continue;
                    else
                    {
                        double cosTheta = latticeBasis.col(0).dot(q2.cartesian()) / (latticeBasis.col(0).norm() * q2.cartesian().norm());
                        if (cosTheta - 1 > -epsilon) cosTheta = 1.0;
                        if (cosTheta + 1 < epsilon) cosTheta = -1.0;
                        double theta = acos(cosTheta);
                        Eigen::Rotation2D<double> rotation(theta);
                        IntScalarType key= theta*keyScale;
                        temp.insert(std::pair<IntScalarType,MatrixDimD>(key,rotation.toRotationMatrix()));
                    }
                }
                else
                {
                    RationalApproximations <IntScalarType> alphaSequence(ratio, maxDen, ratio*maxStrain);
                    for (const auto& alpha: alphaSequence.approximations)
                    {
                        RationalLatticeDirection<dm> q2ByAlpha(Rational<IntScalarType>(alpha.d, alpha.n), q2);

                        mn.col(0) = q2 * alpha.d;
                        md.col(0).setConstant(alpha.n);

                        double s1 =
                                (q2ByAlpha.cartesian().squaredNorm() - latticeBasis.col(0).squaredNorm()) /
                                latticeBasis.col(0).squaredNorm();
                        if (abs(s1) > maxStrain) continue;

                        for (const auto &pair2: coPrimePairs) {
                            VectorDimI qIndices;
                            qIndices << pair2.first, pair2.second;
                            LatticeVector<dm> r2(qIndices, *this);
                            double ratio2= r2.cartesian().norm() / latticeBasis.col(1).norm();
                            RationalApproximations<IntScalarType> betaSequence(ratio2, maxDen,maxStrain*ratio2);
                            for(const auto& beta : betaSequence.approximations) {
                                RationalLatticeDirection<dm> r2ByBeta(Rational<IntScalarType>(beta.d, beta.n), r2);
                                double s2 = (r2ByBeta.cartesian().squaredNorm() - latticeBasis.col(1).squaredNorm()) /
                                            latticeBasis.col(1).squaredNorm();
                                double s3 = (q2ByAlpha.cartesian().dot(r2ByBeta.cartesian()) -
                                             latticeBasis.col(0).dot(latticeBasis.col(1))) /
                                            (latticeBasis.col(0).norm() * latticeBasis.col(1).norm());
                                if (abs(s2) > maxStrain || abs(s3) > maxStrain) continue;

                                // calculate the deformation gradient
                                MatrixDimD F = q2ByAlpha.cartesian() * reciprocalBasis.col(0).transpose() +
                                               r2ByBeta.cartesian() * reciprocalBasis.col(1).transpose();
                                if (F.determinant() < 0) continue;

                                output.push_back( F );
                            }
                        }
                    }
                }
            }
            // sort the angles if rotations are being asked
            if (abs(maxStrain) < epsilon)
                std::transform(temp.begin(), temp.end(),
                               std::back_inserter(output),
                               [](const std::pair<IntScalarType,MatrixDimD>& p) {
                                   return p.second;
                               });
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
 * -# The cross product of two lattice vectors is a reciprocal lattice direction.
 * Similarly, the cross product of two reciprocal lattice vectors is a lattice direction.
 * The cross product is enabled only for dim=3
 *  @snippet testLattice.cpp Cross product
 *
 *
 * Full code:
 */
}
#endif
