//
// Created by Nikhil Chandra Admal on 5/17/24.
//

#ifndef OILAB_GBMESOSTATES_H
#define OILAB_GBMESOSTATES_H

#include <GbShifts.h>
#include <deque>
#include <GbMesoState.h>
#include <GbContinuum.h>

namespace gbLAB {
    /*! Class template that aids in the construction of an ensemble of GB mesostates.
     *
     */
    template<int dim>
    class GbMesoStateEnsemble : public GbShifts<dim>
    {
        using VectorDimD = typename LatticeCore<dim>::VectorDimD;
        using BicrystalLatticeVectors= std::vector<LatticeVector<dim>>;
        using Constraints= std::deque<std::pair<LatticeVector<dim>,VectorDimD>>;

        static std::deque<Constraints> enumerateConstraints(const GbShifts<dim>& gbs);

        /*!
         * \brief Constructs \p bicrystalConfig and \p ensembleGbCslVectors
         * @param gbs - a const reference to GBShifts<dim> object
         * @param ensembleGbCslVectors (output) - the CSL vectors of the ensemble's grain boundary
         * @param scales - an integer array (of size=\p dim) representing the scaling of the mesostate ensemble
         * @return Lattice vectors of lattices \f$\mathcal A\f$ and \f$\mathcal B\f$ in the ensemble's bicrystal
         */
        static BicrystalLatticeVectors getBicrystalConfig(const GbShifts<dim>& gbs,
                                                          std::vector<LatticeVector<dim>>& ensembleGbCslVectors,
                                                          const Eigen::Vector<int,dim>& scales);
    public:
        /*!
         * CSL vectors that define the ensemble's grain boundary region
         */
        std::vector<LatticeVector<dim>> ensembleGbCslVectors;

        /*!
         * A vector of lattice vectors in the ensemble's bicrystal.
         */
        BicrystalLatticeVectors bicrystalConfig;


        /*!
         * A collection of mesostates that form the ensemble
         */
        std::deque<Constraints> constraintsEnsemble;

        GbMesoStateEnsemble(const Gb<dim>& gb,
                            const ReciprocalLatticeVector<dim>& axis,
                            const double& bhalfMax= 0.3,
                            const Eigen::Vector<int,dim>& scales= Eigen::Vector<int,dim>::Ones());

        /*!
         * \brief Constructs an ensemble of mesostates
         * @param gbs - a const reference to GBShifts<dim> object
         * @param scaledGbCslVectors - the CSL vectors of the ensemble's grain boundary
         * @param bicrystalConfig - Lattice vectors of lattices \f$\mathcal A\f$ and \f$\mathcal B\f$
         * @return A deque of mesostates
         */
        std::deque<GbMesoState<dim>> collectMesoStates(const std::string& filename="");

    };

}

#endif //OILAB_GBMESOSTATES_H