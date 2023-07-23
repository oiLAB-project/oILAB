#include <Eigen/Core>
#include <Spectra/GenEigsSolver.h>
#include <iostream>
#include <LatticeModule.h>
#include <Diff.h>

using namespace Spectra;
using namespace gbLAB;

int main()
{
    const int dim=3;
    //const Eigen::array<Eigen::Index,2> n{32,32};
    const Eigen::array<Eigen::Index,3> n{5,5,5};
    Eigen::Matrix<double,dim,dim> A;
    //A << 1.0, 0.0,
    //        0.0, 1.0;
    A << 1.0, 0.0, 0.0,
         0.0, 0.5, 0.0,
         0.0, 0.0, 0.25;
    A=2*M_PI*A;
    Lattice<dim> L(A);
    Diff<dim> dx2({2,0,0},A,n);
    Diff<dim> dy2({0,2,0},A,n);
    Diff<dim> dz2({0,0,2},A,n);
    auto op= dx2+dy2+dz2;
    auto oop= 1*op + (-2)*op + 2*op;

    GenEigsSolver<decltype(oop)> eigs(oop, 20, 45);

    // Initialize and compute
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

