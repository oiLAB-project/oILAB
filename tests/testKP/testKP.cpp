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
#include "KP.h"

using namespace Spectra;
using namespace gbLAB;

/*
 * Example: The 1D Kronig--Penney model
 *
 */

int main()
{
    const int dim=1;
    const int nx= 60;
    Eigen::array<Eigen::Index,dim> n{nx};
    // Lattice system
    Eigen::Matrix<double,dim,dim> A;

    double a= 1.0;
    A << a;
    Lattice<dim> L(A);
    Eigen::Matrix<double,dim,1> basisAtoms;
    Eigen::Vector<double,dim> offset;
    offset << 0;
    basisAtoms.col(0)= offset;

    /*
    double a= 3.0;
    A << a;
    Lattice<dim> L(A);
    Eigen::Matrix<double,dim,3> basisAtoms;
    Eigen::Vector<double,dim> offset;
    offset << 0;
    basisAtoms.col(0)= offset;
    offset << 1;
    basisAtoms.col(1)= offset;
    offset << 2;
    basisAtoms.col(2)= offset;
     */

    MultiLattice<dim> M(A,basisAtoms);

    Diff<dim> dx1({1},A,n);
    Diff<dim> dx2({2},A,n);
    auto laplacian= dx2;

    KP<dim> kp(a);

    // Generate potential V
    PeriodicFunction<double,dim> V= PeriodicFunction<double,dim>::kernelConvolution(n,M,kp);
    Multiplication<dim> mop_V(V);


    std::vector<Eigen::Vector<double,dim>> etaPoints;
    int numberOfKPoints=100;
    for(int i = 0; i<numberOfKPoints; ++i)
        etaPoints.emplace_back((double)i/(a*numberOfKPoints)-1/(2.0*a));


    std::map<int,const Eigen::VectorXd> bands;
    Eigen::MatrixXcd waveFunctionReal, waveFunctionImag;
    waveFunctionReal.resize(n[0],numberOfKPoints);
    waveFunctionImag.resize(n[0],numberOfKPoints);
    std::cout << "Starting band structure calculation" << std::endl;

    ofstream outputBands("bandStructure.txt");
    ofstream outputWavefunctionReal("waveFunctionReal.txt");
    ofstream outputWavefunctionImag("waveFunctionImag.txt");
    ofstream outputPhase("phase.txt");
    #pragma omp parallel for
    for (int i=0; i < etaPoints.size(); i++) {
        auto schrodingerReal = -0.5 * laplacian + mop_V;
        auto schrodingerComplex = 2.0 * M_PI * (etaPoints[i](0) * dx1);
        ComplexOperator<decltype(schrodingerReal), decltype(schrodingerComplex),dim> schrodinger(schrodingerReal,
                                                                                             schrodingerComplex);

        GenEigsSolver<decltype(schrodinger)> eigs(schrodinger, 40, 100);

        eigs.init();
        auto nconv = eigs.compute(SortRule::SmallestReal);

        // Retrieve results
        Eigen::VectorXd evalues;
        Eigen::VectorXcd evector;
        if (eigs.info() == CompInfo::Successful)
        {
            evalues = eigs.eigenvalues().real();

            std::vector<int> indices(evalues.size());
            std::iota(indices.begin(), indices.end(), 0);
            std::sort(indices.begin(), indices.end(),
                      [&](int A, int B) -> bool {
                          return evalues(A) < evalues(B);
                      });

            std::sort(evalues.data(),evalues.data()+evalues.size());
            /*
             * eigenvectors: (v,-iv); (v*,iv*). Sum = (v+v*, -iv+iv*) = 2(Re v, Im v);
             * if v= a+ib,  v* = a,-ib and v+v* = 2 Re v
             * -iv = b-ia)  and iv* = b+ia, -iv+iv* = 2b = 2 Im v
             */
            evector= (eigs.eigenvectors().col(indices[0]) + eigs.eigenvectors().col(indices[1]))/2.0;

            // e = e + 2 * pi^2 * eta^2
            evalues= evalues + Eigen::VectorXd::Constant( evalues.size(), 2 * std::pow(M_PI, 2) * etaPoints[i].squaredNorm());
            std::cout << "\n Number of eigenvalues converged: " << nconv  << std::endl;
        }

        #pragma omp critical
        {
            bands.insert({i,evalues});
            waveFunctionReal.col(i)= evector(Eigen::seqN(0,n[0]));
            waveFunctionImag.col(i)= evector.tail(n[0]);
        }
    }

    const Eigen::IOFormat fmt(15, 0, "\t", "", "\t", "", "", "");
    outputBands << "# Column 1 represents k points, and column i+1 represents the i-th band" << std::endl;
    for (int i=0; i < etaPoints.size(); i++) {
        outputBands << etaPoints[i].transpose().format(fmt) << bands[i].transpose().format(fmt) << std::endl;
    }
    const Eigen::IOFormat fmt2(15, 0, "\t", "\n", "\t", "", "", "");

    outputWavefunctionReal << "# Column i represents the real part of phi_{1 k_i}(x)" << std::endl;
    outputWavefunctionReal << waveFunctionReal.real().format(fmt2) << std::endl;
    outputWavefunctionImag << "# Column i represents the imaginary part of phi_{1 k_i}(x)" << std::endl;
    outputWavefunctionImag << waveFunctionImag.real().format(fmt2) << std::endl;
    outputPhase << (waveFunctionImag.real().array()/waveFunctionReal.real().array()).atan().format(fmt2) << std::endl;

    outputBands.close();
    outputPhase.close();
    outputWavefunctionReal.close();
    outputWavefunctionImag.close();
    return 0;
}

