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
    using IntScalarType = typename LatticeCore<dim>::IntScalarType;
    private:
        MatrixDimI getBasisT(const BiCrystal<dim>& bc, const ReciprocalLatticeDirection<dim>& n);

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

        const Lattice<dim> T;

        const MatrixDimI basisT;

        ReciprocalLatticeVector<dim> getReciprocalLatticeVectorInT(const ReciprocalLatticeVector<dim>& v) const;

        ReciprocalLatticeDirection<dim> getReciprocalLatticeDirectionInT(const ReciprocalLatticeVector<dim>& v) const;



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

        template<int dm=dim>
        typename std::enable_if<dm==2,LatticeVector<dim>>::type
        getPeriodVector(const ReciprocalLatticeVector<dim>& axis) const;

        template<int dm=dim>
        typename std::enable_if<dm==3,LatticeVector<dim>>::type
        getPeriodVector(const ReciprocalLatticeVector<dm>& axis) const;

        /*! This function outputs/prints a grain boundary (two lattices that form the GB,
         * CSL, and the DSCL) bounded by a box defined using
         * input box vectors. The box vectors have to be linearly independent lattice
         * vectors. The \p boxVectors[0] is not parallel to the boundary plane while the
         * remaining box vectors should lie in the GB plane. The function optimizes boxVectors[0]
         * to make the box as orthogonal as possible depending on the \p orthogonality parameter.
         * The length of the box along \p boxVectors[0] is equal to \p 2*boxVectors[0].norm,
         * while the lengths along the remaining box vectors are equal to their norms.
         *
         * The function outputs DSCL lattice points in the GBs neighborhood, which can be controlled
         * by the \p dsclFactor parameter.
         *
         * @tparam dm dimension (int)
         * @param boxVectors three linearly independent lattice vectors. The first box vector
         * is not parallel to the boundary plane, while the remaining box vectors span the GB plane.
         * @param orthogonality (double) a value in the interval \f$[0,1]\f$.
         * @param dsclFactor (int) a factor\f$ \ge 1\f$ to increase the number of outputted DSCL planes
         * @param filename (optional) name of the output file
         * @param orient (optional) While printing to a file, orient the system such that the box sides are along
         * the global x, y, and z axes. The box vectors spanning the grain boundary have to be orthogonal
         * if orient==true. This flag does not
         * influence the returning configuration, only the configuration printed to the file.
         * @return lattice points of the grain boundary bounded by the box (std::vector<LatticeVector<dim>>).
         */
        template<int dm=dim>
        typename std::enable_if<dm==2 || dm==3,std::vector<LatticeVector<dim>>>::type
        box(std::vector<LatticeVector<dim>>& boxVectors,
            const double& orthogonality,
            const int& dsclFactor,
            std::string filename= "",
            bool orient=false) const;

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
 * This examples demonstrates the construction of a grain boundary in a 3D bicrystal.
 * We will construct a [111] tilt \f$\Sigma 273\f$ - STGB (5 6 -11) and an
 * ATGB (17  23 -40)(145  136 -281).
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
