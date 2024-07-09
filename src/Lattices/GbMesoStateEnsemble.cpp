//
// Created by Nikhil Chandra Admal on 5/26/24.
//
#include <GbMesoStateEnsemble.h>

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
                                               ensembleCslVectors)),
            constraintsEnsemble(enumerateConstraints( (const GbShifts<dim>&) *this))
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

        std::cout << "Number of mesostates in the ensemble = " << constraintsEnsemble.size() << std::endl;
        std::cout << std::endl;
        std::cout << "------------------------------" << std::endl;
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
    std::deque<GbMesoState<dim>> GbMesoStateEnsemble<dim>::collectMesoStates(const std::string& filename) const
    {
        std::deque<GbMesoState<dim>> mesoStates;

        int count= -1;
        for(const Constraints& constraints : constraintsEnsemble)
        {
            std::deque<std::pair<LatticeVector<dim>,VectorDimD>> bsPairs;
            for(int i=0; i<this->bShiftPairs.size(); ++i) {
                if (constraints(i) == 1) {
                    const auto& b= this->bShiftPairs[i].first;
                    auto shift= this->bShiftPairs[i].second;
                    bsPairs.push_back(std::make_pair(b,shift));
                }
            }

            count++;
            std::cout << "Constructing mesostate " << count << " of " << constraintsEnsemble.size() << std::endl;
            mesoStates.emplace_back(GbMesoState<dim>(this->gb, this->axis, bsPairs, ensembleCslVectors, bicrystalConfig));
            if (!filename.empty())
                mesoStates.back().box(filename + std::to_string(count));
        }
        return mesoStates;
    }

    /*-------------------------------------*/
    template<int dim>
    std::deque<typename GbMesoStateEnsemble<dim>::Constraints> GbMesoStateEnsemble<dim>::enumerateConstraints(const GbShifts<dim>& gbs)
    {
        std::deque<Constraints> constraintsEnsemble;

        auto tuples= XTuplet::generate_tuples(2,gbs.bShiftPairs.size());
        for (auto& tuple : tuples)
            constraintsEnsemble.push_back(tuple);

        return constraintsEnsemble;
    }

    // function - evolveConstraints(MC parameters, filename)
    // initialize constraints
    //      generate a random evolution of the constraints
    //      construct mesostate
    //      accept/reject
    //  output mesostates as the constraints evolve

    template class GbMesoStateEnsemble<3>;
}
