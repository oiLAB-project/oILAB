//
// Created by Nikhil Chandra Admal on 8/14/24.
//

#ifndef OILAB_CANONICALTP_H
#define OILAB_CANONICALTP_H
#include <EvolutionAlgorithm.h>
#include <utility>
#include <map>

namespace gbLAB {
    // CanonicalTP is an evolution algorithm with a transition probability
    template<typename StateType, typename SystemType>
    class CanonicalTP : public EvolutionAlgorithm<StateType, SystemType, CanonicalTP<StateType,SystemType>> {
    private:
        int countTP;
        double currentEnergy, currentDensity;
    public:
        double temperature;
        std::map<StateType, double> stateEnergyMap;

        CanonicalTP(const double temperature);
        double probability(const std::pair<StateType, SystemType>& proposedState,
                           const std::pair<StateType, SystemType>& currentState) ;

    };

}

#endif //OILAB_CANONICALTP_H
