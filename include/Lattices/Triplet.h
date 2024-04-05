//
// Created by Nikhil Chandra Admal on 2/23/24.
//

#ifndef OILAB_TRIPLET_H
#define OILAB_TRIPLET_H
#include<Eigen/Eigen>

namespace gbLAB {
class Triplet : public Eigen::Vector3i {
    public:
        Triplet() = default;

        Triplet(int i, int j, int k) : Eigen::Vector3i(i, j, k) {}

        Triplet(const Eigen::Vector3i &base) : Eigen::Vector3i(base) {}

        // Define < operator to use std::map<Triplet,std::vector<int>>
        bool operator<(const Triplet &rhs) const {
            if (this->operator()(0) < rhs(0)) return true;
            if (rhs(0) < this->operator()(0)) return false;
            if (this->operator()(1) < rhs(1)) return true;
            if (rhs(1) < this->operator()(1)) return false;
            if (this->operator()(2) < rhs(2)) return true;
            return false;
        }

        virtual ~Triplet() {}
    };

}

#endif //OILAB_TRIPLET_H
