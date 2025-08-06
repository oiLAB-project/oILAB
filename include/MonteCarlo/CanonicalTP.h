//
// Created by Nikhil Chandra Admal on 8/14/24.
//

#ifndef OILAB_CANONICALTP_H
#define OILAB_CANONICALTP_H
#include <EvolutionAlgorithm.h>
#include <utility>
#include <map>
#include <fstream>

namespace gbLAB {
    // CanonicalTP is an evolution algorithm with a transition probability
    template<typename StateType, typename SystemType>
    class CanonicalTP : public EvolutionAlgorithm<StateType, SystemType, CanonicalTP<StateType,SystemType>> {
    private:
        int countTP;
        double currentEnergy, currentDensity;
        std::ofstream output;
        std::string lmpLocation;
        std::string potentialName;
    public:
        double temperature;
        std::map<StateType, double> stateEnergyMap;

        CanonicalTP(const std::string& lmpLocation,
                    const std::string& potentialName,
                    const double& temperature,
                    const std::string& filename="");
        double probability(const std::pair<StateType, SystemType>& proposedState,
                           const std::pair<StateType, SystemType>& currentState) ;

    };

}

#include <CanonicalTPImplementation.h>
#endif //OILAB_CANONICALTP_H
