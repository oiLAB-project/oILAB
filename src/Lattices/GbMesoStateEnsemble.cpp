//
// Created by Nikhil Chandra Admal on 5/26/24.
//
#include <GbMesoStateEnsemble.h>

namespace gbLAB {
    template<int dim>
    GbMesoStateEnsemble<dim>::GbMesoStateEnsemble(const Gb<dim>& gb, const ReciprocalLatticeVector<dim>& axis):
            GbShifts<dim>(gb,axis),
            mesoStates(enumerateMesoStates(gb,axis,this->bShiftPairs))
    {
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
    std::deque<GbMesoState<dim>> GbMesoStateEnsemble<dim>::enumerateMesoStates(const Gb<dim>& gb,
                                                            const ReciprocalLatticeVector<dim>& axis,
                                                            const std::vector<std::pair<LatticeVector<dim>,VectorDimD>>& bShiftPairs)
    {
        std::deque<GbMesoState<dim>> mesoStates;
        const int n= bShiftPairs.size();

        // form k combinations
        for(int k=1; k<=n; ++k) {
            for(const auto& combination : IntegerMath<int>::comb(n,k)) {
                std::deque<std::pair<LatticeVector<dim>, VectorDimD>> bs;
                for (const int& index: combination) {
                    bs.push_back(bShiftPairs[index]);
                }
                mesoStates.emplace_back(GbMesoState<dim>(gb,axis,bs,GbShifts<dim>::getGbCslVectors(gb,axis)));
            }
        }
        return mesoStates;
    }


    //template class GbMesoStateEnsemble<2>;
    template class GbMesoStateEnsemble<3>;
}
