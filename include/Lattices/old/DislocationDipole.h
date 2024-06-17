//
// Created by Nikhil Chandra Admal on 2/17/24.
//

#ifndef OILAB_DISLOCATIONDIPOLE_H
#define OILAB_DISLOCATIONDIPOLE_H
#include "Eigen/Dense"

namespace gbLAB {
    class DislocationDipole
    {
        using Matrix2d= Eigen::Matrix2d;
        using Vector2d= Eigen::Vector2d;

        template <typename T> int sgn(T val) const
        {
            return (T(0) < val) - (val < T(0));
        }
    public:
        //const Matrix2d ends;
        //const Vector2d b;
        //const Vector2d shift;
        //const int nImages;
        Matrix2d ends;
        Vector2d b;
        Vector2d shift;
        int nImages;

        DislocationDipole(const Matrix2d& end_in,const Vector2d& b_in,
                          const Vector2d& shift_in,  const int& nImages_in);

        double solidAngle(const Vector2d& x, const int& branch) const;
        Vector2d solidAngleGradient(const Vector2d& x, const int& branch) const;
        Vector2d displacement(const Vector2d& x, const int& branch) const;
        static std::pair<Vector2d,double> localPosition(const Vector2d& A,const Vector2d& B,const Vector2d& x);

    };
}
#endif //OILAB_DISLOCATIONDIPOLE_H
