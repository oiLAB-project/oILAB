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

        /*! This function generates deformations \f$\mathbf F\f$ such that the deformations of *this lattice share moire supercells
         *  with the undeformed *this lattice
         */
        template<int dm=dim>
        typename std::enable_if<dm==2,std::vector<MatrixDimD>>::type
        generateCoincidentLattices(const double& maxStrain,
                                   const int& maxDen=50,
                                   const int& N=30) const
        {
            std::vector<MatrixDimD> output(generateCoincidentLattices(*this,maxStrain,maxDen,N));
            return output;
        }
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
                                   const int& N=30) const
        {
            std::vector<MatrixDimD> output;
            std::map<IntScalarType,MatrixDimD> temp;
            MatrixDimI mn(MatrixDimI::Zero());
            MatrixDimI md(MatrixDimI::Ones());
            std::vector<std::pair<int,int>> coPrimePairs= farey(N,false);
            double epsilon=1e-8;
            IntScalarType keyScale= 1e6;

            for(const auto& pair1 : coPrimePairs)
            {
                VectorDimI pIndices;
                pIndices << pair1.first, pair1.second;
                LatticeVector<dm> q2(pIndices,undeformedLattice);
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
                            LatticeVector<dm> r2(qIndices, undeformedLattice);
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
        box(std::vector<LatticeVector<dim>> boxVectors, const std::string& filename= "") const
        {
            for(const LatticeVector<dim>& boxVector : boxVectors)
            {
                assert(this == &boxVector.lattice && "Box vectors belong to different lattice.");
            }

            // Form the box lattice
            MatrixDimD C;
            for (int i=0; i<dim; ++i) {
                C.col(i) = boxVectors[i].cartesian();
            }
            assert(C.determinant() != 0);
            Lattice<dim> boxLattice(C);

            auto planeParallelBasis= planeParallelLatticeBasis(boxVectors[1].cross(boxVectors[2]),true);
            MatrixDimD A;
            MatrixDimI mat;
            for(int i=0; i<dim; ++i) {
                A.col(i) = planeParallelBasis[i].cartesian();
                mat.col(i)=planeParallelBasis[i].latticeVector();
            }

            // l - lattice with plane parallel basis, where the second and third basis vectors
            //     lie on the plane spanned by boxVectors[1] and boxVectors[2]
            Lattice<dim> l(A);

            // Form plane parallel vectors v1 and v2
            // v1 is along boxVector[1]
            LatticeDirection<dim> v1(MatrixDimIExt<IntScalarType,dim>::adjoint(mat)*boxVectors[1],l);
            assert(v1.latticeVector()(0)==0);
            IntScalarType x,y;
            IntegerMath<IntScalarType>::extended_gcd(v1.latticeVector()(1),-v1.latticeVector()(2),x,y);
            LatticeVector<dim> temp(l);
            temp << 0,y,x;
            LatticeDirection<dim> v2(temp);

            int areaRatio= round((boxLattice.latticeBasis.col(1).cross(boxLattice.latticeBasis.col(2))).norm()/
                                (l.latticeBasis.col(1).cross(l.latticeBasis.col(2))).norm());
            int scale1= abs(IntegerMath<IntScalarType>::gcd(boxVectors[1]));
            int scale2= round(areaRatio/scale1);
            int scale0= round(abs(C.determinant()/latticeBasis.determinant()/areaRatio));

            std::vector<LatticeVector<dim>> output;

            // r0, r1, and r2 are reciprocal vectors perpendicular to areas spanned by boxVectors[0],
            // boxVectors[1], and boxVectors[2]
            ReciprocalLatticeVector<dim> r0_temp(*this), r1_temp(*this), r2_temp(*this);
            r0_temp= (boxVectors[1].cross(boxVectors[2])).reciprocalLatticeVector();
            if (r0_temp.dot(boxVectors[0]) < 0)
                r0_temp=-1*r0_temp;
            r1_temp= (boxVectors[2].cross(boxVectors[0])).reciprocalLatticeVector();
            if (r1_temp.dot(boxVectors[1]) < 0)
                r1_temp=-1*r1_temp;
            r2_temp= (boxVectors[0].cross(boxVectors[1])).reciprocalLatticeVector();
            if (r2_temp.dot(boxVectors[2]) < 0)
                r2_temp=-1*r2_temp;
            assert(r0_temp.dot(boxVectors[1])==0 && r0_temp.dot(boxVectors[2])==0);
            assert(r1_temp.dot(boxVectors[0])==0 && r1_temp.dot(boxVectors[2])==0);
            assert(r2_temp.dot(boxVectors[0])==0 && r2_temp.dot(boxVectors[1])==0);
            ReciprocalLatticeDirection<dim> r0(r0_temp), r1(r1_temp), r2(r2_temp);

            for(int i=0; i<scale0; ++i)
            {
                for(int j=0; j<scale1; ++j)
                {
                    for(int k=0; k<scale2; ++k) {
                        LatticeVector<dim> vectorIn_l(l);
                        vectorIn_l << i, 0, 0;
                        vectorIn_l= vectorIn_l + j*v1.latticeVector()+k*v2.latticeVector();
                        //LatticeVector<dim> vector(latticeVector(vectorIn_l.cartesian()));
                        LatticeVector<dim> vector((mat*vectorIn_l).eval(),*this);
                        //LatticeDirection<dim> v1(MatrixDimIExt<IntScalarType,dim>::adjoint(mat)*boxVectors[1],l);

                        int c0= IntegerMath<IntScalarType>::positive_modulo(vector.dot(r0),boxVectors[0].dot(r0)) -
                                vector.dot(r0);
                        c0= c0/boxVectors[0].dot(r0);
                        int c1= IntegerMath<IntScalarType>::positive_modulo(vector.dot(r1),boxVectors[1].dot(r1)) -
                                vector.dot(r1);
                        c1= c1/boxVectors[1].dot(r1);
                        int c2= IntegerMath<IntScalarType>::positive_modulo(vector.dot(r2),boxVectors[2].dot(r2)) -
                                vector.dot(r2);
                        c2= c2/boxVectors[2].dot(r2);
                        vector= vector+c0*boxVectors[0]+c1*boxVectors[1]+c2*boxVectors[2];
                        output.push_back(vector);
                    }
                }
            }

            if(!filename.empty()) {
                std::ofstream file;
                file.open(filename);
                if (!file) std::cerr << "Unable to open file";
                file << output.size() << std::endl;
                file << "Lattice=\"";

                for(const auto& vector:boxVectors) {
                    file << vector.cartesian().transpose() << " ";
                }
                file << "\" Properties=atom_types:I:1:pos:R:3" << std::endl;

                for (const auto &vector: output)
                    file << 1 << " " << vector.cartesian().transpose() << std::endl;

                file.close();
            }
            return output;
        }

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
        box(const std::vector<LatticeVector<dim>>& boxVectors, const std::string& filename= "") const
        {
            for(const LatticeVector<dim>& boxVector : boxVectors)
            {
                assert(this == &boxVector.lattice && "Box vectors belong to different lattice.");
            }

            // Form the box lattice
            MatrixDimD C;
            for (int i=0; i<dim; ++i) {
                C.col(i) = boxVectors[i].cartesian();
            }
            assert(C.determinant() != 0);
            Lattice<dim> boxLattice(C);

            // Form plane parallel vectors v0 and v1
            //LatticeDirection<dim> v0(this->latticeDirection(boxVectors[0].cartesian()));
            LatticeDirection<dim> v0(
                    (VectorDimI) boxVectors[0]/IntegerMath<IntScalarType>::gcd(boxVectors[0]),
                    *this);

            ReciprocalLatticeVector<dim> temp(*this);
            temp << v0.latticeVector()(0),-v0.latticeVector()(1);
            LatticeDirection<dim> v1(planeParallelLatticeBasis(temp,true)[0]);

            /*
            IntScalarType x,y;
            IntegerMath<IntScalarType>::extended_gcd(v0.latticeVector()(0),-v0.latticeVector()(1),x,y);
            LatticeVector<dim> temp(*this);
            temp << y,x;

            LatticeDirection<dim> v1(temp);
             */

            int areaRatio= round(abs(boxLattice.latticeBasis.determinant()/
                                            (*this).latticeBasis.determinant()));

            int scale1= abs(IntegerMath<IntScalarType>::gcd(boxVectors[0]));
            int scale2= round(areaRatio/scale1);

            std::vector<LatticeVector<dim>> output;

            // r0 and r1 are reciprocal vectors perpendicular to boxVectors[1] and boxVectors[0], respectively.
            ReciprocalLatticeVector<dim> r1_temp(*this), r0_temp(*this);
            r1_temp << boxVectors[0](1),-boxVectors[0](0);
            if (r1_temp.dot(boxVectors[1]) < 0)
                r1_temp=-1*r1_temp;
            r0_temp << boxVectors[1](1),-boxVectors[1](0);
            if (r0_temp.dot(boxVectors[0]) < 0)
                r0_temp=-1*r0_temp;
            assert(r1_temp.dot(boxVectors[0])==0);
            assert(r0_temp.dot(boxVectors[1])==0);
            ReciprocalLatticeDirection<dim> r1(r1_temp), r0(r0_temp);

            for(int j=0; j<scale1; ++j)
            {
                for(int k=0; k<scale2; ++k) {
                    LatticeVector<dim> vector(j*v0.latticeVector()+k*v1.latticeVector());
                    int c1= IntegerMath<IntScalarType>::positive_modulo(vector.dot(r1),boxVectors[1].dot(r1)) -
                            vector.dot(r1);
                    c1= c1/boxVectors[1].dot(r1);
                    int c2= IntegerMath<IntScalarType>::positive_modulo(vector.dot(r0),boxVectors[0].dot(r0)) -
                            vector.dot(r0);
                    c2= c2/boxVectors[0].dot(r0);
                    vector= vector+c1*boxVectors[1]+c2*boxVectors[0];
                    output.push_back(vector);
                }
            }

            if(!filename.empty()) {
                std::ofstream file;
                file.open(filename);
                if (!file) std::cerr << "Unable to open file";
                file << output.size() << std::endl;
                file << "Lattice=\"";

                for(const auto& vector:boxVectors) {
                    file << vector.cartesian().transpose() << " 0 ";
                }
                file << " 0 0 1 ";
                file << "\" Properties=atom_types:I:1:pos:R:3" << std::endl;

                for (const auto &vector: output)
                    file << 1 << " " << vector.cartesian().transpose() << " 0 " << std::endl;

                file.close();
            }
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
