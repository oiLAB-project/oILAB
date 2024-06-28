//
// Created by Nikhil Chandra Admal on 5/26/24.
//
#include <GbMesoStateEnsemble.h>

namespace gbLAB {
    template<int dim>
    GbMesoStateEnsemble<dim>::GbMesoStateEnsemble(const Gb<dim>& gb,
                                                  const ReciprocalLatticeVector<dim>& axis,
                                                  const double& bhalfMax,
                                                  const Eigen::Vector<int,dim>& scales) :
            GbShifts<dim>(gb,axis,bhalfMax),
            bicrystalConfig(getBicrystalConfig((const GbShifts<dim>&) *this,ensembleGbCslVectors,scales)),
            constraintsEnsemble(enumerateConstraints( (const GbShifts<dim>&) *this))
            //bicrystalConfig(getBicrystalConfig((const GbShifts<dim>&) *this,ensembleGbCslVectors,scales)),
            //mesoStates(collectMesoStates((const GbShifts<dim>&) *this, ensembleGbCslVectors,bicrystalConfig))
    {
        std::cout << "Forming mesostate ensemble with material parameters: " << std::endl;
        std::cout << "lambda = " << GbMaterialTensors::lambda << std::endl;
        std::cout << "mu = " << GbMaterialTensors::mu<< std::endl;

        /*
        for(const auto& mesoState : mesoStates) {
            std::cout << "--------------------------" << std::endl;
            for (const auto& xu: mesoState.xuPairs) {
                std::cout << "x = " << xu.first.transpose() << "; u = " << xu.second.transpose() << std::endl;
            }
        }
        for(const auto& bShift: this->bShiftPairs) {
            std::cout << "--------------------------" << std::endl;
            std::cout << "b = " << bShift.first.cartesian().transpose() << "; s = " << bShift.second.transpose() << std::endl;
            }
        */
    }

    template<int dim>
    typename GbMesoStateEnsemble<dim>::BicrystalLatticeVectors
    GbMesoStateEnsemble<dim>::getBicrystalConfig(const GbShifts<dim>& gbs,
                                                 std::vector<LatticeVector<dim>>& ensembleGbCslVectors,
                                                 const Eigen::Vector<int,dim>& scales)
    {
        std::deque<LatticeVector<dim>> temp;

        for (int i=0; i<dim-1; ++i)
            temp.push_back(scales(i+1)*gbs.gbCslVectors[i]);

        auto basis= gbs.gb.bc.csl.planeParallelLatticeBasis(
                gbs.gb.bc.getReciprocalLatticeDirectionInC(gbs.gb.nA.reciprocalLatticeVector()), true);

        temp.push_front(basis[0].latticeVector());

        while (!temp.empty())
        {
            ensembleGbCslVectors.emplace_back(std::move(temp.front()));
            temp.pop_front();
        }
        gbs.gb.bc.box(ensembleGbCslVectors,1,1);
        ensembleGbCslVectors[0]= scales(0)*ensembleGbCslVectors[0];

        auto allLatticeVectors= gbs.gb.bc.box(ensembleGbCslVectors,1,1);
        std::vector<LatticeVector<dim>> output;

        // include only lattice vectors in A and B
        for (const auto& latticeVector : allLatticeVectors) {
            if (&latticeVector.lattice == &gbs.gb.bc.A || &latticeVector.lattice == &gbs.gb.bc.B)
                output.push_back(latticeVector);
        }

        return output;
    }

    template<int dim>
    std::deque<GbMesoState<dim>> GbMesoStateEnsemble<dim>::collectMesoStates(const std::string& filename)
    {
        // enumerate - form bs pairs
        // split this up into only forming bs. Have a separate static functions to  "collectMesoStates", "evolveMesostates", and "collectMesostates"
        // to be used later for MC simulations
        std::deque<GbMesoState<dim>> mesoStates;

        // form a set of ordered tuple (k1,k2), where ki ranges from 0 to ni-1
        // sz = size of this set
        // form all tuples of dim sz with elements in the range 0 - size(constraintsEnsemble)
        //

        // form k combinations
        // need to work on this. include combinations when the gb is scaled
        int count= -1;
        for(const Constraints& constraints : constraintsEnsemble)
        {
            count++;
            mesoStates.emplace_back(GbMesoState<dim>(this->gb, this->axis, constraints, ensembleGbCslVectors, bicrystalConfig));
            if (!filename.empty())
            {
                std::cout << count << " of " << constraintsEnsemble.size() << std::endl;
                mesoStates.back().box(filename + std::to_string(count));
            }
        }
        return mesoStates;
    }

    template<int dim>
    std::deque<typename GbMesoStateEnsemble<dim>::Constraints> GbMesoStateEnsemble<dim>::enumerateConstraints(const GbShifts<dim>& gbs)
    {
        // enumerate - form bs pairs
        std::deque<Constraints> constraintsEnsemble;
        const int n= gbs.bShiftPairs.size();

        // form k combinations
        for(int k=1; k<=n; ++k) {
            for(const auto& combination : IntegerMath<int>::comb(n,k))
            {
                Constraints bs;
                for (const int& index: combination)
                    bs.push_back(gbs.bShiftPairs[index]);
                constraintsEnsemble.push_back(bs);
            }
        }
        return constraintsEnsemble;
    }
    //template class GbMesoStateEnsemble<2>;
    template class GbMesoStateEnsemble<3>;
}
