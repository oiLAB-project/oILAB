//
// Created by Nikhil Chandra Admal on 12/31/23.
//

#ifndef OILAB_ROTATION_H
#define OILAB_ROTATION_H

#include<Eigen/Eigen>

template<int dim>
class Rotation : public Eigen::Matrix<double,dim,dim>
{
    template<int dm=dim>
    typename std::enable_if<dm==2, Eigen::Matrix<double,dim,dim>  >::type
    static  getMatrix(const Eigen::Matrix<double,dim,dim-1>& orthogonalVectors)
    {
        Eigen::Matrix<double,dim,dim> output;
        output.row(0)= orthogonalVectors.col(0).normalized();
        output.row(1)=Eigen::Rotation2D<double>(M_PI/2)*orthogonalVectors.col(0).normalized();
        return output;
    }
    template<int dm=dim>
    typename std::enable_if<dm==3, Eigen::Matrix<double,dim,dim>  >::type
    static  getMatrix(const Eigen::Matrix<double,dim,dim-1>& orthogonalVectors)
    {
        assert(abs(orthogonalVectors.col(0).dot(orthogonalVectors.col(1))) < FLT_EPSILON);
        Eigen::Matrix<double,dim,dim> output;
        output.row(0)= orthogonalVectors.col(0).normalized();
        output.row(1)= orthogonalVectors.col(1).normalized();
        output.row(2)= (output.row(0).cross(output.row(1))).normalized();
        return output;
    }
public:
    Rotation(const Eigen::Matrix<double,dim,dim-1> orthogonalVectors) :
    Eigen::Matrix<double,dim,dim>(getMatrix(orthogonalVectors))
    {
    }
};
#endif //OILAB_ROTATION_H
