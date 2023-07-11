#include <Eigen/Core>
#include <Spectra/GenEigsSolver.h>
#include <iostream>
#include <FFT.h>
#include <LatticeModule.h>
#include <Laplacian.h>
#include <Multiplication.h>

using namespace Spectra;
using namespace gbLAB;

int main()
{
    Eigen::array<Eigen::Index,2> n{32,32};
    Eigen::Matrix2d A;
    A << 1.0, 0.0,
         0.0, 1.0;
    A=2*M_PI*A;
    Lattice<2> L(A);
    Laplacian<2> op1(A,n);

    PeriodicFunction<double,2> pf(n,L);
    pf.values.setConstant(1);
    Multiplication<double,2> mop(pf);

    auto op= op1+mop;
    GenEigsSolver<decltype(op)> eigs(op, 51, 65);

    eigs.init();
    //int nconv = eigs.compute(SortRule::LargestMagn);
    int nconv = eigs.compute(SortRule::SmallestMagn);

    // Retrieve results
    Eigen::VectorXcd evalues;
    if(eigs.info() == CompInfo::Successful)
        evalues = eigs.eigenvalues();

    std::cout << "Eigenvalues found:\n" << evalues << std::endl;
    std::cout << "\n Number of eigenvalues converged: " << nconv  << std::endl;

    return 0;
}

