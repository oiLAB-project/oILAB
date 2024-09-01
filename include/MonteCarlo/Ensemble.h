//
// Created by Nikhil Chandra Admal on 8/5/24.
//

#ifndef OILAB_ENSEMBLE_H
#define OILAB_ENSEMBLE_H

#include <iostream>

namespace gbLAB {
    template<typename StateType, typename SystemType, typename Derived>
    class Ensemble {
    public:
        const Derived &derived;

        Ensemble() : derived(static_cast<const Derived &>(*this)) {}

        SystemType constructSystem(const StateType &state) const
        {
            return derived.constructSystem(state);
        }

        std::pair<StateType, SystemType> sampleNewState(const StateType &currentState)
        {
            return derived.sampleNewState(currentState);
        }
    };
}

#endif //OILAB_ENSEMBLE_H
