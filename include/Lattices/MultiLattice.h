//
// Created by Nikhil Chandra Admal on 7/4/23.
//

#ifndef OILAB_MULTILATTICE_H
#define OILAB_MULTILATTICE_H

#include<Lattice.h>

namespace gbLAB {
    /*! \brief MultiLattice class
     *
     *  multilattice class description
     * */
    template<int dim>
    class MultiLattice : public Lattice<dim>
    {
    public:
        Eigen::Matrix<double,dim,Eigen::Dynamic> basisAtoms;
        MultiLattice(const Lattice<dim>& _lattice, const Eigen::Matrix<double,dim,Eigen::Dynamic>& _basisAtoms):
            Lattice<dim>(_lattice), basisAtoms(_basisAtoms)
        {}
    };
}
#endif //OILAB_MULTILATTICE_H
