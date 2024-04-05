//
// Created by Nikhil Chandra Admal on 2/18/24.
//

#ifndef OILAB_DISLOCATIONS_H
#define OILAB_DISLOCATIONS_H

#include<DislocationDipole.h>
#include<vector>

namespace gbLAB {
    class Dislocations : private std::vector<DislocationDipole>
    {
        using Matrix2d= Eigen::Matrix2d;
        using Vector2d= Eigen::Vector2d;
        double line_search(const Vector2d& x, const Vector2d& dx, const Vector2d& y, const int& branch) const;
        Matrix2d jacobian(const Vector2d& x, const int& branch) const;
        Vector2d inverseDeformationMap(const Vector2d& x, const int& branch) const;
        const std::vector<DislocationDipole>& dislocations() const;

        public:

        const double a2;
        const Vector2d shift;
        const int nImages;

        Dislocations(const double& a2,
                     const double& shiftSize,
                     const int& nImages);

        void insertDislocationDipole(Matrix2d& ends,Vector2d& b);
        void removeDislocationDipole(const int& index);
        Vector2d deformationMap(const Vector2d& X, const int& branch) const;
        double edgeElasticEnergyKernel(const Vector2d& x1,
                                       const Vector2d& x2,
                                       const Vector2d& b1,
                                       const Vector2d& b2) const;
        double elasticEnergy() const;

    };
}
#endif //OILAB_DISLOCATIONS_H
