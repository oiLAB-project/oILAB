//
// Created by Nikhil Chandra Admal on 12/11/23.
//

#ifndef OILAB_MESOSTATE_H
#define OILAB_MESOSTATE_H

#include <Gb.h>

namespace gbLAB {

    template<int dim>
    class GbShifts
    {
        using VectorDimD = typename LatticeCore<dim>::VectorDimD;
        using VectorDimI = typename LatticeCore<dim>::VectorDimI;
    protected:
        //static std::vector<LatticeVector<dim>> getGbCslVectors(const Gb<dim>& gb, const ReciprocalLatticeVector<dim>& axis);
        static std::vector<std::pair<LatticeVector<dim>,VectorDimD>>  getbShiftPairs(const Gb<dim>& gb,
                                                                                     const std::vector<LatticeVector<dim>>& gbCslVectors,
                                                                                     const double& bhalfMax);

    public:
        const Gb<dim>& gb;
        const ReciprocalLatticeVector<dim>& axis;
        const std::vector<LatticeVector<dim>> gbCslVectors;
        std::vector<std::pair<LatticeVector<dim>,VectorDimD>> bShiftPairs;
        explicit GbShifts(const Gb<dim>& gb,
                          const ReciprocalLatticeVector<dim>& axis,
                          const std::vector<LatticeVector<dim>>& gbCslVectors,
                          const double& bhalfMax= 1);

    };
}
#endif //OILAB_MESOSTATE_H
