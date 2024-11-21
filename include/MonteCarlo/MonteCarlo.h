//
// Created by Nikhil Chandra Admal on 7/29/24.
//

#ifndef OILAB_DOS_H
#define OILAB_DOS_H

#include <vector>
#include <map>
#include <Eigen/Eigen>
#include <randomInteger.h>
#include <cmath>
#include <EvolutionAlgorithm.h>

namespace gbLAB {
    template<typename StateType, typename SystemType, typename EnsembleType, typename EvolveType>
    class MonteCarlo : public EvolutionAlgorithm<StateType, SystemType, EvolveType> {
    public:

        StateType currentState;

        const EnsembleType& ensemble;

        MonteCarlo(const EnsembleType& ensemble, const EvolveType &evolve);

        MonteCarlo(const EnsembleType& ensemble, const EvolveType &evolve, const StateType& state);

        void evolve(const int &maxIterations);
    };


/* ---------------------------------------------------*/


}

#endif //OILAB_DOS_H
