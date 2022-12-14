//
// Created by Nikhil Chandra Admal on 11/5/22.
//

#ifndef OILAB_GB_H
#define OILAB_GB_H

#include "BiCrystal.h"
#include "ReciprocalLatticeDirection.h"

namespace gbLAB
{
    template<int dim>
    class Gb
    {
    using VectorDimI = typename LatticeCore<dim>::VectorDimI ;

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


/*! @example testGenerateGBs.cpp
 * Given a tilt-axis of a 3D lattice, this example demonstrates the construction of multiple tilt GBs of varying
 * misorientations and inclinations
 *
 * -# Define types
 * @snippet testGenerateGBs.cpp Types
 *
 * -# Instantiate a lattice \f$\mathcal A\f$
 * @snippet testGenerateGBs.cpp Lattice
 *
 * -# Specify an axis in the form of a reciprocal lattice direction
 * @snippet testGenerateGBs.cpp Axis
 *
 * -# Generate all rotations \f$\mathbf R\f$ about the given axis that result in a coincidence relation between
 * \f$\mathcal A\f$ and \f$\mathbf R\mathcal A\f$, and form the corresponding bicrystals
 * @snippet testGenerateGBs.cpp Generate bicrystal
 *
 * -# Loop over the generated bicrystals, i.e., loop over misorientation angles
 * @snippet testGenerateGBs.cpp Misorientation
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
 *   -# Output GB properties
 *   @snippet testGenerateGBs.cpp Output
 *
 * Full code:
 */
}

#endif //OILAB_GB_H
