//
// Created by Nikhil Chandra Admal on 7/4/23.
//

#ifndef OILAB_PERIODICFUNCTION_H
#define OILAB_PERIODICFUNCTION_H

#include "Eigen/Dense"
#include "unsupported/Eigen/CXX11/Tensor"
#include <iomanip>

namespace gbLAB {

    template<typename Scalar, int dim>
    class LatticeFunction;

    template<typename Derived, typename Scalar>
    class Function;

    template<typename Scalar, int dim>
    class PeriodicFunction {
    public:
        using dcomplex= std::complex<double>;
        const Eigen::Matrix<double,Eigen::Dynamic,dim> unitCell;
        Eigen::Tensor<Scalar,dim> values;

        // Construct a zero periodic function
        explicit PeriodicFunction(const Eigen::array<Eigen::Index,dim>& n,
                                  const Eigen::Matrix<double,Eigen::Dynamic,dim>&  _unitCell);

        // generates a periodic function from a function centered at the center of a lattice
        template<typename T, typename = T, typename=T, int dm=dim, typename = std::enable_if_t<dm==1>>
        PeriodicFunction(const Eigen::array<Eigen::Index,dim>& n,
                         const Eigen::Matrix<double,Eigen::Dynamic,dim>&  _unitCell,
                         const Function<T,Scalar>& fun);

        template<typename T, typename = T, int dm=dim, typename = std::enable_if_t<dm==2>>
        PeriodicFunction(const Eigen::array<Eigen::Index,dim>& n,
                         const Eigen::Matrix<double,Eigen::Dynamic,dim>&  _unitCell,
                         const Function<T,Scalar>& fun);

        template<typename T, int dm=dim, typename = std::enable_if_t<dm==3>>
        PeriodicFunction(const Eigen::array<Eigen::Index,dim>& n,
                         const Eigen::Matrix<double,Eigen::Dynamic,dim>&  _unitCell,
                         const Function<T,Scalar>& fun);

        LatticeFunction<dcomplex,dim> fft() const;

        double dot(const PeriodicFunction<Scalar,dim>& other) const;

        template <typename T>
        PeriodicFunction<Scalar,dim> kernelConvolution(const Function<T,Scalar>& kernel);
    };

    template<typename Scalar, int dim>
    std::basic_ostream<char>& operator<<(std::basic_ostream<char>& s, const PeriodicFunction<Scalar, dim>& fun)
    {
        s << "x coord, y coord, scalar" << std::endl;
        auto n= fun.values.dimensions();
        for (int i = 0; i < n[0]; i++) {
            for (int j = 0; j < n[1]; j++) {
                Eigen::Vector<double,Eigen::Dynamic> x= i*fun.unitCell.col(0)/n[0] +
                                             j*fun.unitCell.col(1)/n[1];
                const Eigen::IOFormat fmt(15, 0, ", ", "", "\t", "", "", "");
                s << x.transpose().format(fmt) << ", " << std::setw(25) << std::setprecision(15) << fun.values(i,j) << std::endl;
            }
        }
        return s;
    }

}

#endif //OILAB_PERIODICFUNCTION_H
