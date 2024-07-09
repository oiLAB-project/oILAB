//
// Created by Nikhil Chandra Admal on 5/26/24.
//

#ifndef OILAB_GBMESOSTATE_H
#define OILAB_GBMESOSTATE_H

#include <Gb.h>
#include <GbContinuum.h>
#include <LatticeCore.h>
#include <OrderedTuplet.h>

namespace gbLAB {

    /*! Class template that defines a GB mesostate.
     *
     */
    template<int dim>
    class GbMesoState : public GbContinuum<dim> {
        using VectorDimD = typename LatticeCore<dim>::VectorDimD;
        using BicrystalLatticeVectors= std::vector<LatticeVector<dim>>;

        /*!
         * \brief Returns the cartesian coordinates of the CSL vectors that define a mesostate's GB.
         * @param mesoStateCslVectors - a vector (size = \p dim) of CSL vectors that define the box of the mesostate.
         * @return Cartesian coordinates of the \p dim-1 grain boundary CSL vectors
         */
         // ensure that the input is of the right dimension
        static Eigen::Matrix<double, dim,dim-1> getMesoStateGbDomain(const std::vector<LatticeVector<dim>>& mesoStateCslVectors);


        /*!
         * \brief Returns the displacement constraints \f$\textbf u(\textbf x^i)=\textbf u^i\f$ for the GB continuum model from the translation-shift pairs
         * defining the mesostate.
         * @param gb - grain boundary
         * @param bs - a deque of pairs (translation vector \p b, shift vector \p s). \p b is a DSCL vector while \p s is expressed in Cartesian coordinates
         * @return A map between the DSCL coordinates \f$\underline{\textbf x}^i_{\mathcal D}\f$ of \f$\textbf x^i\f$ and Cartesian coordinates of \f$\textbf u^i\f$.
         */
        static std::map<OrderedTuplet<dim+1>, VectorDimD> get_xuPairs(const Gb<dim>& gb,
                                                                      const std::vector<LatticeVector<dim>>& mesoStateCslVectors,
                                                                      const std::deque<std::pair<LatticeVector<dim>, VectorDimD>>& bs);

        /*!
         * \brief Returns the Cartesian coordinates of the lattice vectors of the mesostates's bicrystal in the form of a map that maps the
         * DSCL integer coordinates to the Cartesian coordinates.
         * @param bicrystalConfig - a vector of lattice vectors in the mesostates' bicrystal.
         * @return A map between the DSCL coordinates of lattice vectors in bicrystalConfig to their Cartesian coordinates.
         */
        static std::map<OrderedTuplet<dim+1>,VectorDimD> bicrystalCoordsMap(const Gb<dim>& gb, const BicrystalLatticeVectors& bicrystalConfig);


        /*!
         * \brief Returns the discretization of the grain boundary domain
         * @param mesoStateCslVectors - a vector (size = \p dim) of CSL vectors that define the box of the mesostate.
         * @return an integer array of size \p dim-1
         */
        static std::array<Eigen::Index,dim-1> discretize(const std::vector<LatticeVector<dim>>& mesoStateCslVectors, const Gb<dim>& gb);

    public:

        /*!
         * Grain boundary
         */
        const Gb<dim>& gb;

        /*!
         * Grain boundary tilt axis
         */
        const ReciprocalLatticeVector<dim>& axis;

        /*!
         * A vector (size = \p dim) of CSL vectors that define the box of the mesostate. The second and third vectors
         * should be parallel to the grain boundary, while the third vector should be out of the grain boundary  plane.
         */
         // need the ensure that the inputs respect the above constraint
        const std::vector<LatticeVector<dim>>& mesoStateCslVectors;

        /*!
         * Lattice vectors of lattices \f$\mathcal A\f$ and \f$\mathcal B\f$ that are contained in the box described by mesoStateCslVectors.
         * They are used to construct a function basis for the continuum solution.
         */
        const BicrystalLatticeVectors& bicrystalConfig;


        /*!
         * @param bs a deque of pairs (translation vector \f$\textbf b\f$, shift vector \f$\textbf s\f$) that defines a mesostate. Translating lattice \f$\mathcal A\f$
         *           by \f$\textbf b\f$ results in a CSL shift of \f$\textbf s\f$.
         */
        const std::deque<std::pair<LatticeVector<dim>, VectorDimD>> bs;

         // ensure that b in bs pair belongs to the DSCL vectors and s belongs to the box
         // ensure that the lattice vectors in bicrystalConfig lie entirely in the box
        explicit GbMesoState(const Gb<dim>& gb,
                             const ReciprocalLatticeVector<dim>& axis,
                             const std::deque<std::pair<LatticeVector<dim>,VectorDimD>>& bs,
                             const std::vector<LatticeVector<dim>>& mesoStateCslVectors,
                             const BicrystalLatticeVectors& bicrystalConfig);

        double energy() const;

        /*! This function outputs/prints a grain boundary mesostate
         * @param filename name of the file to be written to
         */
        typename std::enable_if<dim==3,void>::type box(const std::string& filename) const;
    };
}
#endif //OILAB_GBMESOSTATE_H
