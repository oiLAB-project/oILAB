#include <Eigen/Core>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/SymEigsSolver.h>
#include <iostream>
#include <LatticeModule.h>
#include <MultiLattice.h>
#include <Diff.h>
#include <fstream>

using namespace Spectra;
using namespace gbLAB;


int main()
{
    const int dim=3;
    const int nx= 12;
    const int ny= 12;
    const int nz= 12;
    Eigen::array<Eigen::Index,dim> n{nx,ny,nz};

    // Lattice system
    Eigen::Matrix<double,dim,dim> A;
    A << 0.0,   -2.13042249331,   0.0,
         2.46,  -1.23,            0.0,
         0.0,    0.0,             16.0;
    Lattice<dim> L(A);

    // input function: cos(px) sin(qy) cos(rz)
    const int p= 2, q= 3, r= 5;

    const int ndx= 2;
    const int ndy= 3;
    const int ndz= 4;




    Diff<dim> op({ndx,ndy,ndz},A,n);


    // Generate potential V
    PeriodicFunction<dcomplex,dim> V_complex= PeriodicFunction<dcomplex,dim>::kernelConvolution(n,M,mayer);
    PeriodicFunction<double,dim> V(n,M);
    V.values= V_complex.values.real();


    return 0;
}

