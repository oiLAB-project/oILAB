#include <Eigen/Core>
#include <LatticeModule.h>
#include <PeriodicFunction.h>
#include <fstream>
#include "Mayer.h"

using namespace gbLAB;

int main()
{
    const int dim= 3;
    const int nx= 12;
    const int ny= 12;
    const int nz= 12;
    Eigen::array<Eigen::Index,dim> n{nx,ny,nz};
    double tol= 1e-10;
    double maxError= 0;
    // Read the kernel function data

    std::string kernelFilename= "meyer.txt";
    double kernelFunctionDomain= 16.0;
    std::cout << "Domain of the kernel function = " << kernelFunctionDomain << std::endl;

    // Lattice system
    Eigen::Matrix<double,dim,dim> A;
    A << 0.0,   -2.13042249331,   0.0,
         2.46,  -1.23,            0.0,
         0.0,    0.0,             16.0;
    A= A/0.5291772109;
    Lattice<dim> L(A);
    Eigen::Matrix<double,dim,2> basisAtoms;
    Eigen::Vector<double,dim> offset;
    //offset << 0,0;
    offset << 0,0,0;
    basisAtoms.col(0)= offset;
    basisAtoms.col(1)= A.col(0)/3 + 2*A.col(1)/3+offset;
    MultiLattice<dim> M(A,basisAtoms);

    // Generate potential V

    Mayer<dim> mayer(kernelFunctionDomain);

    // Generate potential V
    PeriodicFunction<double,dim> V= PeriodicFunction<double,dim>::kernelConvolution(n,M,mayer);

    // Check for error
    std::cout << "checking for error" << std::endl;
    for (int p=0; p<nz; p++){
        for (int l= 0; l<nx; l++) {
            for (int m = 0; m < ny; m++) {
                Eigen::Vector<double, dim> center = A.col(0) * l / nx +
                                                    A.col(1) * m / ny +
                                                    A.col(2) * p / nz;
                double potential = 0;

                for (int i = -10; i < 10; i++) {
                    for (int j = -10; j < 10; j++) {
                        for (int k = -3; k < 3; k++) {
                            for (const auto &atom: basisAtoms.colwise()) {
                                Eigen::Vector<double, dim> position = i * A.col(0) +
                                                                      j * A.col(1) +
                                                                      k * A.col(2) +
                                                                      atom;
                                //potential += mayer(center - position).real();
                                potential += mayer(center - position);
                            }
                        }
                    }
                }

                double error = abs(V.values(l, m, p) - potential);
                maxError = max(error, maxError);
                try {
                    if (error > tol) {
                        std::cout << "error = " << error << std::endl;
                        throw std::runtime_error("TestPeriodicFunction failed");
                    }
                }
                catch (std::runtime_error &e) {
                    std::cout << e.what() << std::endl;
                    return -1;
                }
            }
        }
    }
    std::cout << "Maximum error = " << maxError << std::endl;


    return 0;
}

