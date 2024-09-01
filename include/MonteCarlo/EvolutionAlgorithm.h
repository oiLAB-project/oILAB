//
// Created by Nikhil Chandra Admal on 8/14/24.
//

#ifndef OILAB_EVOLUTIONALGORITHM_H
#define OILAB_EVOLUTIONALGORITHM_H

#include <utility>

namespace gbLAB {
    template<typename StateType, typename SystemType, typename TransitionProbabilityType>
    class EvolutionAlgorithm {
    public:
        TransitionProbabilityType& transitionProbability;

        EvolutionAlgorithm();

        bool acceptMove(const std::pair<StateType,SystemType>& proposedStateSystem,
                        const std::pair<StateType,SystemType>& currentStateSystem) const;
    };
}
/* ---------------------------------------------------*/
#endif //OILAB_EVOLUTIONALGORITHM_H
