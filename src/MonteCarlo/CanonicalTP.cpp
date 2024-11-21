//
// Created by Nikhil Chandra Admal on 8/14/24.
//
#include <CanonicalTP.h>
#include <algorithm>
#include <OrderedTuplet.h>
#include <GbMesoStateEnsemble.h>

namespace gbLAB {

    template<typename StateType, typename SystemType>
    CanonicalTP<StateType,SystemType>::CanonicalTP(const double& temperature,
                                                   const std::string& filename) : temperature(temperature),
                                                                                  countTP(0)
                                                                                  {
                                                                                      if (!filename.empty())
                                                                                          output.open(filename);
                                                                                  }

    template<typename StateType, typename SystemType>
    double CanonicalTP<StateType,SystemType>::probability(const std::pair<StateType,SystemType>& proposedStateSystem,
                                                          const std::pair<StateType,SystemType>& currentStateSystem)
    {
        // calculate energies here
        const auto& currentState= currentStateSystem.first;
        const auto& currentSystem= currentStateSystem.second;
        const auto& proposedState= proposedStateSystem.first;
        const auto& proposedSystem= proposedStateSystem.second;

        if(countTP==0) {
            const auto &temp = currentSystem.densityEnergy();
            currentDensity = temp.first;
            currentEnergy = temp.second;
            stateEnergyMap[currentState] = currentEnergy;
            //std::cout << "density = " << currentDensity << ", energy = " << currentEnergy << std::endl;
        }
        else {
            assert(stateEnergyMap.find(currentState) != stateEnergyMap.end());
            currentEnergy= stateEnergyMap[currentState];
        }

        double proposedEnergy, proposedDensity;

        //StateType proposedState(proposedStateSystemPair.first);
        //SystemType proposedSystem(proposedStateSystemPair.second);
        if (stateEnergyMap.find(proposedState) != stateEnergyMap.end()) {
            proposedEnergy = stateEnergyMap.at(proposedState);
        }
        else {
            const auto& temp= proposedSystem.densityEnergy();
            proposedDensity= temp.first;
            proposedEnergy= temp.second;
            //proposedEnergy = proposedSystem.energy();
            //std::cout << "density = " << proposedDensity << ", energy = " << proposedEnergy << std::endl;
            stateEnergyMap[proposedState] = proposedEnergy;
        }

        countTP++;
        double delta = proposedEnergy - currentEnergy;

        if(output.is_open())
            output << currentState << "  " << currentState.density() << "     " << currentEnergy << std::endl;
        return std::min(1.0, exp(-delta / temperature));
    }

    template class CanonicalTP<XTuplet,GbMesoState<3>>;
/* ---------------------------------------------------*/
}
