/* This file is part of gbLAB.
 *
 * gbLAB is distributed without any warranty under the MIT License.
 */


#ifndef gbLAB_LatticeCore_cpp_
#define gbLAB_LatticeCore_cpp_

#include <LatticeCore.h>
#include <BestRationalApproximation.h>
#include <Eigen/Dense>

namespace gbLAB
{



template <int dim>
typename LatticeCore<dim>::VectorDimI LatticeCore<dim>::rationalApproximation(VectorDimD nd)
{
    Eigen::Array<IntScalarType, dim, 1> nums = VectorDimI::Zero();

    if (nd.squaredNorm() > 0.0)
    {
        const Eigen::Array<double, dim, 1> nda(nd.array().abs()); // vector of close-to-integer numbers corresponding to lattice coordinates
        size_t maxID = 0;
        const double maxVal(nda.maxCoeff(&maxID));
        nd /= maxVal; // make each value of nd in [-1:1]

        nums = VectorDimI::Ones();
        Eigen::Array<IntScalarType, dim, 1> dens = VectorDimI::Ones();
        IntScalarType denProd = 1;

        for (int k = 0; k < dim; ++k)
        {
            //BestRationalApproximation bra(nd(k), 10000);
            BestRationalApproximation bra(nd(k), 1000);

            nums(k) = bra.num;
            dens(k) = bra.den;
            denProd *= bra.den;
        }

        for (int k = 0; k < dim; ++k)
        {
            nums(k) *= (denProd / dens(k));
        }
    }

    return nums.matrix();
}

template <int dim>
typename LatticeCore<dim>::VectorDimI LatticeCore<dim>::integerCoordinates(const VectorDimD& d,const MatrixDimD& invA)
{
    const VectorDimD nd(invA*d);
    const VectorDimD rd(nd.array().round());
    if ((nd - rd).norm() > roundTol)
    {
        std::cout << "nd=" << nd.transpose() << std::endl;
        std::cout << "rd=" << rd.transpose() << std::endl;
        std::cout << "rounding error = |nd-rd| = " << (nd-rd).norm() << std::endl;
        throw(std::runtime_error("Input vector is not a lattice vector"));
    }
    return rd.template cast<IntScalarType>();
}



    template struct LatticeCore<1>;
    template struct LatticeCore<2>;
    template struct LatticeCore<3>;
    template struct LatticeCore<4>;
    template struct LatticeCore<5>;
} // end namespace
#endif
