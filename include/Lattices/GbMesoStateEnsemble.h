//
// Created by Nikhil Chandra Admal on 5/17/24.
//

#ifndef OILAB_GBMESOSTATES_H
#define OILAB_GBMESOSTATES_H

#include <GbShifts.h>
#include <deque>
#include <GbMesoState.h>
#include <GbContinuum.h>
#include <Ensemble.h>

namespace gbLAB {
    /*! Class template that aids in the construction of an ensemble of GB mesostates.
     *
     */
    template<int dim>
    class GbMesoStateEnsemble : public GbShifts<dim>,
                                public Ensemble<XTuplet,GbMesoState<dim>,GbMesoStateEnsemble<dim>>
    {
        using VectorDimD = typename LatticeCore<dim>::VectorDimD;
        using BicrystalLatticeVectors= std::vector<LatticeVector<dim>>;
        //using Constraints= Eigen::Tensor<int,dim>;
        using Constraints= XTuplet;

        static std::deque<Constraints> enumerateConstraints(const GbShifts<dim>& gbs);

        /*!
         * \brief Constructs \p bicrystalConfig and \p ensembleCslVectors
         * @param gbs - a const reference to GBShifts<dim> object
         * @param ensembleCslVectors (output) - the CSL vectors of the ensemble's grain boundary
         * @param scales - an integer array (of size=\p dim) representing the scaling of the mesostate ensemble
         * @return Lattice vectors of lattices \f$\mathcal A\f$ and \f$\mathcal B\f$ in the ensemble's bicrystal
         */
        static BicrystalLatticeVectors getBicrystalConfig(const GbShifts<dim>& gbs,
                                                          std::vector<LatticeVector<dim>>& ensembleCslVectors);
                                                          //const Eigen::Vector<int,dim>& scales);

        static std::deque<std::tuple<LatticeVector<dim>,VectorDimD,int>> bsPairsFromConstraints(const std::vector<std::pair<LatticeVector<dim>, VectorDimD>>& bShiftPairs,
                                                                                                const Constraints& constraints);
    public:
        /*!
         * CSL vectors that define the ensemble's grain boundary region
         */
        std::vector<LatticeVector<dim>> ensembleCslVectors;

        /*!
         * A vector of lattice vectors in the ensemble's bicrystal.
         */
        BicrystalLatticeVectors bicrystalConfig;


        GbMesoStateEnsemble(const Gb<dim>& gb,
                            const ReciprocalLatticeVector<dim>& axis,
                            std::vector<LatticeVector<dim>>& ensembleCslVectors,
                            const double& bhalfMax);

        /*!
         * \brief Constructs an ensemble of mesostates
         * @param filename-
         * @return A deque of mesostates
         */
        std::deque<GbMesoState<dim>> collectMesoStates(const std::string& filename="") const;


        /*!
         * \brief Evove mesostates using a Monte Carlo algorithm
         * @param filename-
         * @return A deque of mesostates
         */
        std::map<Constraints,GbMesoState<dim>> evolveMesoStates(const double& temperature, const int& resetEvery, const int& maxIterations, const std::string& filename="") const;


        GbMesoState<dim> constructMesoState(const Constraints& constraints) const;

        Constraints sampleNewState(const Constraints& currentConstraints,
                                   const bool& randomize= false) const;

        Constraints initializeState() const;

    };

}

#endif //OILAB_GBMESOSTATES_H