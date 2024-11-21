//
// Created by Nikhil Chandra Admal on 8/6/24.
//
#include <MonteCarlo.h>
#include <OrderedTuplet.h>
#include <GbMesoStateEnsemble.h>
#include <LandauWangTP.h>
#include <CanonicalTP.h>

namespace gbLAB {
    template<typename StateType, typename SystemType, typename EnsembleType, typename EvolveType>
    MonteCarlo<StateType,SystemType,EnsembleType,EvolveType>::
        MonteCarlo(const EnsembleType& ensemble,
                   const EvolveType& evolve) : EvolutionAlgorithm<StateType,SystemType,EvolveType>(evolve),
                                               ensemble(ensemble),
                                               currentState(ensemble.sampleNewState(ensemble.initializeState(), true))
        {}
    template<typename StateType, typename SystemType, typename EnsembleType, typename EvolveType>
    MonteCarlo<StateType,SystemType,EnsembleType,EvolveType>::
    MonteCarlo(const EnsembleType& ensemble,
               const EvolveType& evolve,
               const StateType& state) : EvolutionAlgorithm<StateType,SystemType,EvolveType>(evolve),
                                         ensemble(ensemble),
                                         //currentState(ensemble.sampleNewState(state, false))
                                         currentState(state)
    {}

    template<typename StateType, typename SystemType, typename EnsembleType, typename EvolveType>
    void MonteCarlo<StateType,SystemType,EnsembleType,EvolveType>::evolve(const int& maxIterations)
    {
        int acceptCount = 0;

        for (int i = 0; i < maxIterations; ++i) {
            auto proposedState= ensemble.sampleNewState(currentState, false);
            bool transition = this->acceptMove(std::make_pair(proposedState,ensemble.constructMesoState(proposedState)),
                                               std::make_pair(currentState,ensemble.constructMesoState(currentState)));
            if (transition)
            {
                currentState = proposedState;
                acceptCount++;
            }
        }
    }

    template class MonteCarlo<XTuplet,GbMesoState<3>,GbMesoStateEnsemble<3>,LandauWangTP<XTuplet,GbMesoState<3>>>;
    template class MonteCarlo<XTuplet,GbMesoState<3>,GbMesoStateEnsemble<3>,CanonicalTP<XTuplet,GbMesoState<3>>>;

}
