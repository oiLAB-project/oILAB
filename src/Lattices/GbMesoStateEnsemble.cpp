//
// Created by Nikhil Chandra Admal on 5/26/24.
//
#define PY_SSIZE_T_CLEAN
#include <GbMesoStateEnsemble.h>
#include <randomInteger.h>

namespace gbLAB {
    template<int dim>
    GbMesoStateEnsemble<dim>::GbMesoStateEnsemble(const Gb<dim>& gb,
                                                  const ReciprocalLatticeVector<dim>& axis,
                                                  //const Eigen::Vector<int,dim>& scales) :
                                                  std::vector<LatticeVector<dim>>& ensembleCslVectors,
                                                  const double& bhalfMax):
            GbShifts<dim>(gb,
                          axis,
                          std::vector<LatticeVector<dim>>(ensembleCslVectors.begin()+1,
                                                          ensembleCslVectors.end()),
                          bhalfMax),
            ensembleCslVectors(ensembleCslVectors),
            bicrystalConfig(getBicrystalConfig((const GbShifts<dim>&) *this,
                                               ensembleCslVectors))
    {
        std::cout << "--------------------GBMesoStateEnsemble class construction ---------------------------" << std::endl;
        std::cout << "Forming mesostate ensemble with material parameters: ";
        std::cout << "lambda = " << GbMaterialTensors::lambda;
        std::cout << "; mu = " << GbMaterialTensors::mu<< std::endl;
        std::cout << std::endl;

        std::cout << "Ensemble CSL vectors:"  << std::endl;
        for(const auto& latticeVector : ensembleCslVectors)
            std::cout << latticeVector.cartesian().transpose() <<  std::endl;
        std::cout << std::endl;

    }

    /*-------------------------------------*/
    template<int dim>
    typename GbMesoStateEnsemble<dim>::BicrystalLatticeVectors
    GbMesoStateEnsemble<dim>::getBicrystalConfig(const GbShifts<dim>& gbs,
                                                 std::vector<LatticeVector<dim>>& ensembleCslVectors)
    {
        auto allLatticeVectors= gbs.gb.bc.box(ensembleCslVectors,1,1,"bc.txt");
        std::vector<LatticeVector<dim>> output;

        // include only lattice vectors in A and B
        for (const auto& latticeVector : allLatticeVectors) {
            if (&latticeVector.lattice == &gbs.gb.bc.A || &latticeVector.lattice == &gbs.gb.bc.B)
                output.push_back(latticeVector);
        }

        return output;
    }
    /*-------------------------------------*/
    template<int dim>
    std::deque<std::tuple<LatticeVector<dim>,typename GbMesoStateEnsemble<dim>::VectorDimD,int>>
        GbMesoStateEnsemble<dim>::bsPairsFromConstraints(const std::vector<std::pair<LatticeVector<dim>,typename GbMesoStateEnsemble<dim>::VectorDimD>>& bShiftPairs,
                                                                                                                        const typename GbMesoStateEnsemble<dim>::Constraints& constraints)
    {
        using VectorDimD = typename GbMesoStateEnsemble<dim>::VectorDimD;
        assert(bShiftPairs.size() == constraints.size());

        std::deque<std::tuple<LatticeVector<dim>,VectorDimD,int>> bsPairs;
        for(int i=0; i<bShiftPairs.size(); ++i) {
            if (constraints(i) == 1 || constraints(i) == 2) {
                const auto& b= bShiftPairs[i].first;
                auto shift= bShiftPairs[i].second;
                bsPairs.push_back(std::make_tuple(b,shift,constraints(i)));
            }
        }
        return bsPairs;
    }

    /*-------------------------------------*/
    template<int dim>
    std::map<typename GbMesoStateEnsemble<dim>::Constraints,GbMesoState<dim>> GbMesoStateEnsemble<dim>::collectMesoStates(const std::string& filename) const
    {
        std::deque<Constraints> constraintsEnsemble(enumerateConstraints( (const GbShifts<dim>&) *this));
        std::cout << "Number of mesostates in the ensemble = " << constraintsEnsemble.size() << std::endl;
        std::cout << std::endl;
        std::cout << "------------------------------" << std::endl;
        std::map<Constraints,GbMesoState<dim>> mesoStates;

        int count= -1;
        for(const Constraints& constraints : constraintsEnsemble)
        {
            std::deque<std::tuple<LatticeVector<dim>,VectorDimD,int>> bsPairs(bsPairsFromConstraints(this->bShiftPairs,constraints));
            try {
                //mesoStates.emplace_back(constructMesoState(constraints));
                mesoStates.emplace(constraints,constructMesoState(constraints));
                count++;
                std::cout << "Constructing mesostate " << count << " of " << constraintsEnsemble.size() << std::endl;
                std::cout << "Mesostate signature:  " << constraints.transpose() << std::endl;
                if (!filename.empty())
                    //mesoStates.back().box(filename + std::to_string(count));
                    mesoStates.at(constraints).box(filename + std::to_string(count));
            }
            catch(std::runtime_error& e)
            {
                std::cout << e.what() << std::endl;
            }

        }
        return mesoStates;
    }

    /*-------------------------------------*/
    template<int dim>
    GbMesoState<dim> GbMesoStateEnsemble<dim>::constructMesoState(const Constraints& constraints) const {
        std::deque<std::tuple<LatticeVector<dim>,VectorDimD,int>> bsPairs(bsPairsFromConstraints(this->bShiftPairs,constraints));
        try {
            GbMesoState<dim> mesostate(this->gb, this->axis, bsPairs, ensembleCslVectors, bicrystalConfig);
            return mesostate;
        }
        catch(std::runtime_error& e)
        {
            throw(e);
        }
    }

    /*-------------------------------------*/
    template<int dim>
    std::map<typename GbMesoStateEnsemble<dim>::Constraints, GbMesoState<dim>> GbMesoStateEnsemble<dim>::evolveMesoStates(const double& temperature, const int& resetEvery, const int& maxIterations, const std::string& filename) const
    {
        std::map<Constraints,GbMesoState<dim>> mesoStates;
        std::map<Constraints,double> visitedStatesEnergy;

        // create the initial random mesostate
        Constraints initialConstraints(this->bShiftPairs.size());
        initialConstraints.setZero();
        std::cout << "Initial mesostate signature = " << initialConstraints.transpose() << std::endl;
        GbMesoState<dim> initial_ms(this->gb,
                                    this->axis,
                                    bsPairsFromConstraints(this->bShiftPairs,initialConstraints),
                                    ensembleCslVectors,
                                    bicrystalConfig);
        //double newEnergy= initial_ms.energy();
        std::pair<double, double> newDensityEnergyPair= initial_ms.densityEnergy();
        double newEnergy= newDensityEnergyPair.second;
        int count= 0;
        int mesoStateCount= 0;
        double currentEnergy= 0.0;
        Constraints currentConstraints(initialConstraints);
        bool transition= true;

        for(int i=0; i< maxIterations; ++i) {
            // form mesostate with currentConstraints
            if (transition) {
                GbMesoState<dim> current_ms(this->gb,
                                            this->axis,
                                            bsPairsFromConstraints(this->bShiftPairs, currentConstraints),
                                            ensembleCslVectors,
                                            bicrystalConfig);
                currentEnergy = newEnergy;

                // box the mesostate only if it is a new state
                if (mesoStates.find(currentConstraints) == mesoStates.end()) {
                    mesoStates.insert({currentConstraints, current_ms});
                    if (!filename.empty())
                        current_ms.box(filename + std::to_string(mesoStateCount));
                    std::cout << "Step " << i << " Accept " << count << " Mesostate " << mesoStateCount << " constraints = "
                              << currentConstraints.transpose() << "; energy = " << currentEnergy << std::endl;
                    mesoStateCount++;
                }
                transition = false;
            }

            // new mesostate construction
            bool randomize= false;
            if (count % resetEvery == 0 && count != 0) {
                std::cout << "Starting over using ramdomized constraints" << std::endl;
                randomize= true;
            }
            /*
            Constraints newConstraints(this->bShiftPairs.size());
            bool msConstructionSuccess= false;
            while(!msConstructionSuccess) {
                try {
                    newConstraints= currentConstraints;
                    // alter the constraints
                    if (count % resetEvery == 0 && count != 0) {
                        for (auto &elem: newConstraints)
                            elem = random<int>(0, 2);
                    }
                   else {
                        int randomSpot = random<int>(0, this->bShiftPairs.size() - 1);
                        newConstraints(randomSpot) = random<int>(0, 2);
                    }
                    GbMesoState<dim> temp(this->gb,
                                          this->axis,
                                          bsPairsFromConstraints(this->bShiftPairs, newConstraints),
                                          ensembleCslVectors,
                                          bicrystalConfig);
                   msConstructionSuccess= true;
                }
                catch(std::runtime_error& e)
                {
                    //throw(e);
                }
            }
            GbMesoState<dim> new_ms(this->gb,
                                    this->axis,
                                    bsPairsFromConstraints(this->bShiftPairs, newConstraints),
                                    ensembleCslVectors,
                                    bicrystalConfig);
                                    */
            /*
            auto newConstraintsMesoStatePair= sampleNewMesoState(currentConstraints,randomize);
            Constraints newConstraints(newConstraintsMesoStatePair.first);
            GbMesoState<dim> new_ms(newConstraintsMesoStatePair.second);
             */
            Constraints newConstraints(sampleNewState(currentConstraints,randomize));
            GbMesoState<dim> new_ms(constructMesoState(newConstraints));

            if (visitedStatesEnergy.find(newConstraints) != visitedStatesEnergy.end())
                newEnergy = visitedStatesEnergy.at(newConstraints);
            else {
                //newEnergy = new_ms.energy();
                newDensityEnergyPair= new_ms.densityEnergy();
                newEnergy= newDensityEnergyPair.second;
                visitedStatesEnergy[newConstraints] = newEnergy;
            }

            double delta = newEnergy - currentEnergy;
            double probability = std::min(0.9, exp(-delta / temperature));

            if (random<double>(0.0, 1.0) <= probability) {
                currentConstraints = newConstraints;
                count++;
                transition = true;
            }
        }
        return mesoStates;
    }
    /*-------------------------------------*/
    template<int dim>
    std::deque<typename GbMesoStateEnsemble<dim>::Constraints> GbMesoStateEnsemble<dim>::enumerateConstraints(const GbShifts<dim>& gbs)
    {
        std::deque<Constraints> constraintsEnsemble;

        auto tuples= XTuplet::generate_tuples(3,gbs.bShiftPairs.size());
        for (auto& tuple : tuples) {
            if (!(tuple.array()==0).all())
                constraintsEnsemble.push_back(tuple);
        }

        return constraintsEnsemble;
    }

    /*-------------------------------------*/
    template<int dim>
    typename GbMesoStateEnsemble<dim>::Constraints GbMesoStateEnsemble<dim>::sampleNewState(const Constraints& currentConstraints,
                                                                                            const bool& randomize) const
    {
        // new mesostate construction
        Constraints newConstraints(this->bShiftPairs.size());
        bool msConstructionSuccess = false;
        while (!msConstructionSuccess) {
            try {
                newConstraints = currentConstraints;
                // alter the constraints
                if (randomize){
                    for (auto &elem: newConstraints)
                        elem = random<int>(0, 2);
                } else {
                    int randomSpot = random<int>(0, this->bShiftPairs.size() - 1);
                    newConstraints(randomSpot) = random<int>(0, 2);
                }
                GbMesoState<dim> temp(this->gb,
                                      this->axis,
                                      bsPairsFromConstraints(this->bShiftPairs, newConstraints),
                                      ensembleCslVectors,
                                      bicrystalConfig);
                msConstructionSuccess = true;
            }
            catch (std::runtime_error &e) {
                //throw(e);
            }
        }
        return newConstraints;
    }

        /*-------------------------------------*/
    template<int dim>
    typename GbMesoStateEnsemble<dim>::Constraints GbMesoStateEnsemble<dim>::initializeState() const
    {
        Constraints initialConstraints(this->bShiftPairs.size());
        initialConstraints.setZero();
        return initialConstraints;
    }


    template class GbMesoStateEnsemble<3>;
}
