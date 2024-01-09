//
// Created by Nikhil Chandra Admal on 11/5/22.
//

#ifndef OILAB_GB_H
#define OILAB_GB_H

#include "BiCrystal.h"
#include "ReciprocalLatticeDirection.h"
#include "Rotation.h"

namespace gbLAB
{
    template<int dim>
    class Gb
    {
    using VectorDimI = typename LatticeCore<dim>::VectorDimI;
    using VectorDimD = typename LatticeCore<dim>::VectorDimD;
    using MatrixDimD = typename LatticeCore<dim>::MatrixDimD;
    using MatrixDimI = typename LatticeCore<dim>::MatrixDimI;

    public:
        /*!
         * Bicrystal formed by two lattices, say \f$\mathcal A\f$ and \f$\mathcal B\f$
         */
        const BiCrystal<dim>& bc;
        /*!
         * GB normal described with respect to lattice \f$\mathcal A\f$
         */
        const ReciprocalLatticeDirection<dim> nA;
        /*!
         * GB normal described with respect to lattice \f$\mathcal B\f$
         */
        const ReciprocalLatticeDirection<dim> nB;

        /*!
         * \brief Computes the step height of a disconnection formed by displacing lattice \f$\mathcal A\f$
         * by a Burgers vector \f$\textbf d\f$.
         * @param d - Burgers vector that belongs to the DSCL of the bicrystal
         * @return step height
         */
        double stepHeightA(const LatticeVector<dim>& d) const;
        /*!
         * \brief Computes the step height of a disconnection formed by displacing lattice \f$\mathcal B\f$
         * by a Burgers vector \f$\textbf d\f$.
         * @param d - Burgers vector that belongs to the DSCL of the bicrystal
         * @return step height
         */
        double stepHeightB(const LatticeVector<dim>& d) const;

        /*!
         * \brief Constructs a grain boundary of a given orientation in a bicrystal
         * @param bc - Bicrystal formed by two lattices \f$\mathcal A\f$ and \f$\mathcal B\f$.
         * @param n - Reciprocal lattice direction in dual lattices \f$\mathcal A^*\f$ or \f$\mathcal B^*\f$.
         */
        Gb(const BiCrystal<dim>& bc, const ReciprocalLatticeDirection<dim>& n);

        /*! This function outputs/prints a grain boundary (two lattices that form the GB,
         * CSL, and the DSCL) bounded by a box defined using
         * input box vectors. The box vectors have to be linearly independent lattice
         * vectors. The @param boxVectors[0] is not parallel to the boundary plane while the
         * remaining box vectors lie on the GB plane. The function optimizes boxVectors[0]
         * to make the box as orthogonal as possible depending on the @param orthogonality parameter.
         * The length of the box along boxVectors[0] is equal to 2*boxVectors[0].norm,
         * while the lengths along the remaining box vectors are equal to their norms.
         *
         * The function outputs only those DSCL lattice points in a neighborhood of the GB.
         * The parameter @param dsclFactor can be used to increase the number of outputted DSCL planes.
         *
         * @tparam dm dimension (int)
         * @param boxVectors three linearly independent lattice vectors. The first box vector
         * is not parallel to the boundary plane, while the remaining box vectors span the GB plane.
         * @param orthogonality (double) a value in the interval \f$[0,1]\f$.
         * @param dsclFactor (int) a factor>1 to increase the number of outputted DSCL planes
         * @param filename (optional) name of the output file
         * @param orient (optional) While printing to a file, orient the system such that the box sides are along
         * the global x, y, and z axes. The box vectors have to be orthogonal if orient==true. This flag does not
         * influence the returning configuration, only the configuration printed to the file.
         * @return lattice points of the grain boundary bounded by the box (std::vector<LatticeVector<dim>>).
         */
        template<int dm=dim>
        typename std::enable_if<dm==2 || dm==3,std::vector<LatticeVector<dim>>>::type
        box(std::vector<LatticeVector<dim>> boxVectors,
            const double& orthogonality,
            const int& dsclFactor,
            std::string filename= "",
            bool orient=false) const
        {
            assert(orthogonality>=0.0 && orthogonality<=1.0 &&
            "The \"orthogonality\" parameter should be between 0.0 and 1.0");
            assert(dsclFactor>=1 &&
            "The \"dsclFactor\" should be greater than 1.");
            assert(boxVectors.size()==dim);
            for(const auto& boxVector : boxVectors)
            {
                assert(&bc.csl == &boxVector.lattice &&
                        "Box vectors do not belong to the CSL.");
            }
            for (auto iter= std::next(boxVectors.begin()); iter < boxVectors.end(); iter++)
                assert((*iter).dot(bc.getReciprocalLatticeDirectionInC(nA.reciprocalLatticeVector())) == 0 &&
                       "Box vectors not parallel to the grain boundary.");

            // Structure matrix of the inputted box vectors
            MatrixDimD C;
            for (int i=0; i<dim; ++i) {
                C.col(i) = boxVectors[i].cartesian();
            }
            assert(abs(C.determinant()) > FLT_EPSILON && "Box volume is equal to zero.");


            // Adjust boxVector[0] such that it is as orthogonal as possible to the GB
            auto boxVectorTemp= boxVectors[0];
            auto nC= bc.getReciprocalLatticeDirectionInC(nA.reciprocalLatticeVector());
            std::cout << "Input box height= " << abs(nC.planeSpacing()*boxVectors[0].dot(nC)) << std::endl;
            auto basis= bc.csl.planeParallelLatticeBasis(nC,true);
            int planesToExplore= nC.stacking();
            std::cout << "Exploring " << planesToExplore << " planes" << std::endl;
            MatrixDimI boxLatticeIndices;
            boxLatticeIndices.col(0)= boxVectors[0];
            for (int i=1; i<dim; ++i)
                boxLatticeIndices.col(i)= boxVectors[i]/IntegerMath<long long int>::gcd(boxVectors[i]);
            double minDotProduct= M_PI/2;
            int minStep;

            auto boxVectorUpdated(boxVectors[0]);

            for(int i=0;i<planesToExplore;++i)
            {
                boxVectorTemp= boxVectors[0]+i*basis[0].latticeVector();
                boxLatticeIndices.col(0)= boxVectorTemp;
                Lattice<dim> boxLattice(bc.csl.latticeBasis*boxLatticeIndices.template cast<double>());
                ReciprocalLatticeVector<dim> rC(boxLattice);
                rC(0)=1;
                VectorDimI temp=
                        boxLatticeIndices*boxLattice.planeParallelLatticeBasis(ReciprocalLatticeDirection<dim>(rC),true)[0].latticeVector();
                boxVectorTemp= LatticeVector<dim>(temp,bc.csl);
                double dotProduct= abs(acos(boxVectorTemp.cartesian().normalized().dot(rC.cartesian().normalized())));
                if(dotProduct<minDotProduct) {
                    minDotProduct = dotProduct;
                    boxVectorUpdated= boxVectorTemp;
                    if (dotProduct < (1-orthogonality)*M_PI/2)
                        break;
                }

            }
            boxVectors[0]=boxVectorUpdated;
            C.col(0)= boxVectors[0].cartesian();
            std::cout << "Updated box height= " << abs(nC.planeSpacing()*boxVectors[0].dot(nC)) << std::endl;


            // form the rotation matrix used to orient the system
            MatrixDimD rotation= Eigen::Matrix<double,dim,dim>::Identity();;
            Eigen::Matrix<double,dim,dim-1> orthogonalVectors;
            if (orient) {
                if (dim==3) {
                    orthogonalVectors.col(0) = C.col(1).normalized();
                    orthogonalVectors.col(1) = C.col(2).normalized();
                }
                else if (dim==2)
                    orthogonalVectors.col(0) = C.col(1).normalized();

                rotation=Rotation<dim>(orthogonalVectors);
            }
            assert((rotation*rotation.transpose()).template isApprox(Eigen::Matrix<double,dim,dim>::Identity())
                    && "Cannot orient the grain boundary. Box vectors are not orthogonal.");


            std::vector<LatticeVector<dim>> configurationA, configurationB, configurationC, configurationD;
            std::vector<LatticeVector<dim>> configuration;

            std::vector<LatticeVector<dim>> boxVectorsInA, boxVectorsInB, boxVectorsForCSL(boxVectors), boxVectorsInDSCL;
            boxVectorsForCSL[0]=2*boxVectorsForCSL[0];

            for(const auto& boxVector : boxVectors) {
                boxVectorsInA.push_back(bc.A.latticeVector(boxVector.cartesian()));
                boxVectorsInB.push_back(bc.B.latticeVector(boxVector.cartesian()));
                auto dsclVector=bc.getLatticeDirectionInD(boxVector).latticeVector();
                dsclVector= round(boxVector.cartesian().norm()/dsclVector.cartesian().norm())*dsclVector;
                boxVectorsInDSCL.push_back(dsclVector);
            }
            auto dsclVector=bc.getLatticeDirectionInD(boxVectors[0]).latticeVector();
            auto nD= bc.getReciprocalLatticeDirectionInD(nC.reciprocalLatticeVector());
            dsclVector= dsclFactor*dsclVector;
            if(abs((dsclFactor*dsclVector).dot(nD)) < abs(boxVectorsInDSCL[0].dot(nD)))
                boxVectorsInDSCL[0]= dsclFactor*dsclVector;


            std::vector<LatticeVector<dim>> boxVectorsForDSCL(boxVectorsInDSCL);
            boxVectorsForDSCL[0]=2*boxVectorsInDSCL[0];
            auto nboxVectorsInA= boxVectorsInA;
            auto nboxVectorsInB= boxVectorsInB;
            nboxVectorsInA[0]=-1*boxVectorsInA[0];
            nboxVectorsInB[0]=-1*boxVectorsInB[0];

            if (boxVectors[0].dot(bc.getReciprocalLatticeDirectionInC(nA.reciprocalLatticeVector())) > 0)
            {
                configurationB= bc.B.box(boxVectorsInB);
                configurationA= bc.A.box(nboxVectorsInA);
            }
            else
            {
                configurationA= bc.A.box(boxVectorsInA);
                configurationB= bc.B.box(nboxVectorsInB);
            }

            LatticeVector<dim> origin(-1*boxVectors[0]);
            configurationC= bc.csl.box(boxVectorsForCSL);
            configurationD= bc.dscl.box(boxVectorsForDSCL);
            for(auto& vector : configurationC)
                vector= vector + origin;
            for(auto& vector : configurationD)
                vector= vector + LatticeVector<dim>(-1*boxVectorsInDSCL[0]);

            configuration= configurationA;
            configuration.insert(configuration.end(),configurationB.begin(),configurationB.end());
            configuration.insert(configuration.end(),configurationC.begin(),configurationC.end());
            configuration.insert(configuration.end(),configurationD.begin(),configurationD.end());
            std::cout << "Number of lattice points in A    = " << configurationA.size() << std::endl;
            std::cout << "Number of lattice points in B    = " << configurationB.size() << std::endl;
            std::cout << "Number of lattice points in CSL  = " << configurationC.size() << std::endl;
            std::cout << "Number of lattice points in DSCL = " << configurationD.size() << std::endl;

            if(!filename.empty()) {
                std::ofstream file;
                file.open(filename);
                if (!file) std::cerr << "Unable to open file";
                file << configurationA.size() + configurationB.size() + configurationC.size() + configurationD.size() << std::endl;
                file << "Lattice=\"";


                if (dim==2) {
                    file << (rotation*boxVectorsForCSL[0].cartesian()).transpose() << " " << 0.0 << "  ";
                    file << (rotation*boxVectorsForCSL[1].cartesian()).transpose() << " " << 0.0 << "  ";
                    file << " 0 0 1";
                    file << "\" Properties=atom_types:I:1:pos:R:3:radius:R:1 PBC=\"F T T\" origin=\"";
                    file << (rotation*origin.cartesian()).transpose() << " 0.0\"" << std::endl;
                    for (const auto &vector: configurationA)
                        file << 1 << " " << (rotation*vector.cartesian()).transpose() << " " << 0.0 << "  " << 0.05 << std::endl;
                    for (const auto &vector: configurationB)
                        file << 2 << " " << (rotation*vector.cartesian()).transpose() << " " << 0.0 << "  " << 0.05 << std::endl;
                    for (const auto &vector: configurationC)
                        file << 3 << " " << (rotation*vector.cartesian()).transpose() << " " << 0.0 << "  " << 0.2 << std::endl;
                    for (const auto &vector: configurationD)
                        file << 4 << " " << (rotation*vector.cartesian()).transpose() << " " << 0.0 << "  " << 0.01 << std::endl;
                }
                else if (dim==3){
                    file << (rotation*boxVectorsForCSL[0].cartesian()).transpose()  << " ";
                    file << (rotation*boxVectorsForCSL[1].cartesian()).transpose() << " ";
                    file << (rotation*boxVectorsForCSL[2].cartesian()).transpose();
                    file << "\" Properties=atom_types:I:1:pos:R:3:radius:R:1 PBC=\"F T T\" origin=\"";
                    file << (rotation*origin.cartesian()).transpose() << "\"" << std::endl;

                    for (const auto &vector: configurationA)
                        file << 1 << " " << (rotation*vector.cartesian()).transpose() << "  " << 0.05 << std::endl;
                    for (const auto &vector: configurationB)
                        file << 2 << " " << (rotation*vector.cartesian()).transpose() << "  " << 0.05 << std::endl;
                    for (const auto &vector: configurationC)
                        file << 3 << " " << (rotation*vector.cartesian()).transpose() << "  " << 0.2 << std::endl;
                    for (const auto &vector: configurationD)
                        file << 4 << " " << (rotation*vector.cartesian()).transpose() << "  " << 0.01 << std::endl;
                }


                file.close();
            }
            return configuration;
        }

    };

/*!
 * @example testGb.cpp
 * This example demonstrates the computation of step heights and the use of Gb class
 *
 * -# Initialize two lattices
 *  @snippet testGb.cpp Lattice
 *
 * -# Form the bicrystal
 *  @snippet testGb.cpp Bicrystal
 *
 * -# Define the grain boundary normal (w.r.t lattice 2) and form the GB
 *  @snippet testGb.cpp GB
 *
 * -# Define the Burgers vector \f$b\f$ of a disconnection and compute the two step heights,
 * \f$h_{\mathcal A}\f$ and \f$h_{\mathcal B}\f$.
 *  @snippet testGb.cpp Step Height
 *
 * -# Compute the CSL plane spacing \f$H\f$ parallel to the grain boundary
 *  @snippet testGb.cpp Plane Spacing
 *
 * -# Check the conditions: \f$ \mod(h_{\mathcal A} - h_{\mathcal B} - \textbf{b} \cdot \hat{\textbf{n}}_{\mathcal A}, H) = 0 \f$, where
 * \f$\hat{\textbf{n}}_{\mathcal A}\f$ and \f$\hat{\textbf{n}}_{\mathcal B}\f$ are the outward unit normals to lattices \f$\mathcal A\f$ and \f$\mathcal B\f$, respectively.
 *  @snippet testGb.cpp Check
 *
 *
 * Full code:
*/

/*!
 * @example testGb3d.cpp
 * This example demonstrates the computation of step heights and the use of Gb class
 *
 * -# Initialize the first lattice L1
 *  @snippet testGb3d.cpp lattice1
 *
 * -# Define the tilt axis \f$[1\,1\,1]\f$ as a reciprocal vector in L1
 *  @snippet testGb3d.cpp axis
 *
 * -# Form the second lattice L2 by rotating L1 about the tilt axis
 *  @snippet testGb3d.cpp lattice2
 *
 * -# Define the grain boundary normal (w.r.t lattice 1) and form the GB
 *  @snippet testGb3d.cpp gb
 *
 * -# To print, construct a grain boundary bounded by a box formed by three CSL vectors -
 * period vector and vectors along the tilt axis and the grain boundary normal. Period vector
 * is the shortest CSL vector normal to the tilt axis and the grain boundary normal.
 *
 *  -# Calculate the period vector (periodC).
 *  glideA is a vector in lattice L1, and periodC is the shortest CSL vector along glideA.
 *  @snippet testGb3d.cpp period vector
 *
 *  -# Calculate the CSL vector along the normal (normalC).
 *  normalA is a vector in lattice L1, and normalC is the shortest CSL vector along normalA.
 *  @snippet testGb3d.cpp normal vector
 *
 *  -# Calculate the CSL vector along the tilt axis (vectorAlongAxisC).
 *  vectorAlongAxisA is a vector in lattice L1, and vectorAlongAxisC is the shortest CSL vector along vectorAlongAxisA.
 *  @snippet testGb3d.cpp axis vector
 *
 * -# Form the box vectors and output the grain boundary to a file
 *  @snippet testGb3d.cpp box vectors
 * Full code:
*/

/*! @example testGenerateGBs.cpp
 * Given a tilt-axis of a 3D lattice, this example demonstrates the construction of multiple tilt GBs of varying
 * misorientations and inclinations, and the calculation of grain boundaries' disconnection modes
 *
 * -# Define types
 * @snippet testGenerateGBs.cpp Types
 *
 * -# Instantiate a lattice \f$\mathcal A\f$
 * @snippet testGenerateGBs.cpp Lattice
 *
 * -# Specify the Cartesian coordinates of an axis. This should be parallel to a reciprocal lattice vector of
 * \f$\mathcal A\f$
 * @snippet testGenerateGBs.cpp Axis
 *
 * -# Generate all rotations \f$\mathbf R\f$ about the given axis that result in a coincidence relation between
 * \f$\mathcal A\f$ and \f$\mathbf R\mathcal A\f$, and form the corresponding bicrystals
 * @snippet testGenerateGBs.cpp Generate bicrystal
 *
 * -# Loop over the generated bicrystals, i.e., loop over misorientation angles
 * @snippet testGenerateGBs.cpp Misorientation
 *
 * -# SNF output for each misorientation
 * @snippet testGenerateGBs.cpp SNF
 *
 * -# Output the invariance property of the CSL in the following steps:
 * 1) compute reduced basis vectors \f$\mathbf d_1\f$ and \f$\mathbf d_2\f$ for the DSCL,
 * 2) compute the CSL shifts \f$\mathbf s_1\f$ and \f$\mathbf s_2\f$ when lattice \f$\mathcal A\f$ is
 * displaced by \f$\mathbf d_1\f$ and \f$\mathbf d_2\f$, respectively.
 * @snippet testGenerateGBs.cpp Invariance
 *
 *  -# Generate grain boundaries of varying inclinations
 *  @snippet testGenerateGBs.cpp  Generate GBs
 *
 *  -# Loop over inclination angles
 *  @snippet testGenerateGBs.cpp Inclination
 *
 *   -# For each GB, compute its period and the Burgers vector of the glide disconnection
 *   @snippet testGenerateGBs.cpp Glide
 *
 *   -# Note the reference GB with respect to which inclination angles are measured
 *   @snippet testGenerateGBs.cpp Reference
 *
 *   -# Output properties of GBs whose period is less than 100 Angstrom
 *   @snippet testGenerateGBs.cpp Output
 *
 *
 *
 * Full code:
 */
}

#endif //OILAB_GB_H
