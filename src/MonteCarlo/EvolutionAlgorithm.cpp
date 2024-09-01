//
// Created by Nikhil Chandra Admal on 8/14/24.
//
#include <EvolutionAlgorithm.h>
#include <randomInteger.h>
#include <OrderedTuplet.h>
#include <CanonicalTP.h>
#include <LandauWangTP.h>
#include <GbMesoState.h>

namespace gbLAB {
    template<typename StateType, typename SystemType, typename TransitionProbabilityType>
    EvolutionAlgorithm<StateType, SystemType, TransitionProbabilityType>::EvolutionAlgorithm() :
        transitionProbability(static_cast<TransitionProbabilityType &>(*this))
    { }

    template<typename StateType, typename SystemType, typename TransitionProbabilityType>
    bool EvolutionAlgorithm<StateType, SystemType, TransitionProbabilityType>::acceptMove(const std::pair<StateType,SystemType>& proposedStateSystem,
                                                                                          const std::pair<StateType,SystemType>& currentStateSystem) const
    {
        double probability = transitionProbability.probability(proposedStateSystem,currentStateSystem);
        if (random<double>(0.0, 1.0) <= probability)
            return true;
        else
            return false;
    }

    template class EvolutionAlgorithm<XTuplet,GbMesoState<3>,CanonicalTP<XTuplet,GbMesoState<3>>>;
    template class EvolutionAlgorithm<XTuplet,GbMesoState<3>,LandauWangTP<XTuplet,GbMesoState<3>>>;
}
