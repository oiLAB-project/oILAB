#include <Eigen/Core>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/SymEigsSolver.h>
#include <iostream>
#include <FFT.h>
#include <LatticeModule.h>
#include <Diff.h>
#include <Multiplication.h>
#include <fstream>
#include <ComplexOperator.h>
#include "Mayer.h"

using namespace Spectra;
using namespace gbLAB;


int main()
{
    const int dim=3;
    const int nev= 30;
    const int ncv= 100;
    const int nx= 36;
    const int ny= 36;
    const int nz= 36;
    const int div= 8;
    Eigen::array<Eigen::Index,dim> n{nx,ny,nz};
    // Lattice system
    Eigen::Matrix<double,dim,dim> A;
    double interlayerSpacing= 3.335;
    A << 0.0,   -2.13042249331,   0.0,
         2.46,  -1.23,            0.0,
         0.0,    0.0,             2*interlayerSpacing;
    A= A/0.5291772109;
    Lattice<dim> L(A);
    Eigen::Matrix<double,dim,4> basisAtoms;

    double interlayerSpacingBohr= interlayerSpacing/0.5291772109;
    basisAtoms.col(0).setZero();
    //basisAtoms(2,0)= -interlayerSpacingBohr/2.0;
    basisAtoms(2,0)= 0.0;
    basisAtoms.col(1)= A.col(0)/3 + 2*A.col(1)/3;
    //basisAtoms(2,1)= -interlayerSpacingBohr/2.0;
    basisAtoms(2,1)= 0.0;

    basisAtoms.col(2).setZero();
    //basisAtoms(2,2)= interlayerSpacingBohr/2.0;
    basisAtoms(2,2)= interlayerSpacingBohr;

    Eigen::Matrix<double,dim,dim> R;
    //R= Eigen::AngleAxis<double>(M_PI/3,Eigen::Vector3d::UnitZ());
    R= Eigen::AngleAxis<double>(0.0,Eigen::Vector3d::UnitZ());
    basisAtoms.col(3)= R*basisAtoms.col(1);
    //basisAtoms(2,3)= interlayerSpacingBohr/2.0;
    basisAtoms(2,3)= interlayerSpacingBohr;

    MultiLattice<dim> M(A,basisAtoms);
    Diff<dim> dx1({1,0,0},A,n);
    Diff<dim> dy1({0,1,0},A,n);
    Diff<dim> dz1({0,0,1},A,n);
    Diff<dim> dx2({2,0,0},A,n);
    Diff<dim> dy2({0,2,0},A,n);
    Diff<dim> dz2({0,0,2},A,n);
    auto laplacian= dx2+dy2+dz2;


    Mayer<dim> mayer(16.0);

    // Generate potential V
    PeriodicFunction<double,dim> V= PeriodicFunction<double,dim>::kernelConvolution(n,M,mayer);
    ofstream potential("potential.csv");
    potential << V << std::endl;
    potential.close();
    Multiplication<dim> mop_V(V);


    // make this part of the lattice class
    // PointM is along a reciprocal direction
    Eigen::Vector<double,dim> pointM= 0.5*L.reciprocalBasis.col(0);
    // PointK is along a lattice direction
    Eigen::Vector<double,dim> pointK= L.latticeBasis.col(0).normalized()*pointM.norm()*2/sqrt(3);
    std::cout << "point K: " << pointK.transpose() << std::endl;
    std::cout << "point M: " << pointM.transpose() << std::endl;
    std::vector<Eigen::Vector<double,dim>> etaPoints;
    std::cout << "number of divisions = " << div << std::endl;
    // G - K - M - G
    for (int j=0; j < div; j++){
        etaPoints.emplace_back(j * pointK / div);
    }
    for (int j=0; j<div; j++){
        etaPoints.emplace_back(pointK+j*(pointM-pointK)/div);
    }
    for (int j=0; j<div; j++){
        etaPoints.emplace_back(pointM-j*pointM/div);
    }

    std::map<int,const Eigen::VectorXd> bands;
    std::cout << "Starting band structure calculation" << std::endl;

    ofstream output("bandStructure.txt");
    #pragma omp parallel for
    for (int i=0; i < etaPoints.size(); i++) {
        auto schrodingerReal = -0.5 * laplacian + mop_V;
        auto schrodingerComplex = 2.0 * M_PI * (etaPoints[i](0) * dx1 +
                                                                      etaPoints[i](1) * dy1 +
                                                                      etaPoints[i](2) * dz1);
        ComplexOperator<decltype(schrodingerReal), decltype(schrodingerComplex),dim> schrodinger(schrodingerReal,
                                                                                             schrodingerComplex);

        GenEigsSolver<decltype(schrodinger)> eigs(schrodinger, nev, ncv);

        eigs.init();
        auto nconv = eigs.compute(SortRule::SmallestReal);

        // Retrieve results
        Eigen::VectorXd evalues;
        if (eigs.info() == CompInfo::Successful) {
            evalues = eigs.eigenvalues().real();
        std::sort(evalues.data(),evalues.data()+evalues.size());
        evalues= evalues + Eigen::VectorXd::Constant( evalues.size(), 2 * std::pow(M_PI, 2) * etaPoints[i].squaredNorm());
        std::cout << "\n Number of eigenvalues converged: " << nconv  << std::endl;
        }

        #pragma omp critical
        {
            bands.insert({i,evalues});
        }
    }

    for (int i=0; i < etaPoints.size(); i++) {
        const Eigen::IOFormat fmt(15, 0, "\t", "", "\t", "", "", "");
        output << i << etaPoints[i].transpose().format(fmt) << bands[i].transpose().format(fmt) << std::endl;
    }
    output.close();
    return 0;
}

