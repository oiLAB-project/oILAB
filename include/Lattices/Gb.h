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
        const BiCrystal<dim>& bc;
        const ReciprocalLatticeDirection<dim> nA;
        const ReciprocalLatticeDirection<dim> nB;
//        const double cslPlaneSpacing;
//        const double dsclPlaneSpacing;

        double stepHeightA(const LatticeVector<dim>& d) const;
        double stepHeightB(const LatticeVector<dim>& d) const;

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
}

#endif //OILAB_GB_H
