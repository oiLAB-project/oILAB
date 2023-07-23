#include <Eigen/Core>
#include <Spectra/GenEigsSolver.h>
#include <iostream>
#include <LatticeModule.h>
#include <Operator.h>
#include <ComplexOperator.h>
#include <fstream>

using namespace Spectra;
using namespace gbLAB;


namespace gbLAB {
    template<int dim>
    class Mat : public Operator<Mat<dim>,dim> {

    public:
        using Derived= Operator<Mat<dim>,dim>;
        const Eigen::MatrixXd& M;

        explicit Mat(const Eigen::MatrixXd& M_, const Eigen::array<Eigen::Index,dim>& n_) :
                        Derived(Eigen::MatrixXd::Zero(dim,dim),n_),
                        M(M_) { }

        void perform_op(const double *x_in, double *y_out) const {
            Eigen::Map<const Eigen::VectorXd> x(x_in, this->rows());
            Eigen::Map<Eigen::VectorXd> y(y_out, this->rows());
            y = M * x;
        }
    };
}


int main()
{
    const int dim= 2;
    const int nx= 10, ny= 10;

    const int matSize= nx*ny;
    const int k= 40;
    const double tol= 1e-7;

    Eigen::MatrixXd BMat(matSize,matSize);
    Eigen::MatrixXd CMat(matSize,matSize);
    Eigen::VectorXd evaluesOfA(k);

    std::string Bdata   = "realMatrix.txt";
    std::string Cdata   = "complexMatrix.txt";
    std::string evaluesData= "evalues.txt";
    std::cout << "Reading the operators " << std::endl;
    std::ifstream file1(Bdata);
    std::ifstream file2(Cdata);
    std::ifstream file3(evaluesData);

    if(!file1 || !file2 || !file3)
    {
        // Print an error and exit
        std::cerr << "ERROR: files could not be opened for reading!" << std::endl;
        exit(1);
    }


    for (int i=0; i<k; i++)
        file3 >> evaluesOfA(i);


    for(int j=0; j<matSize; j++)
        for(int i=0; i<matSize; i++){
            file1 >> BMat(i,j);
            file2 >> CMat(i,j);
        }

    Mat<2> B(BMat,{nx,ny});
    Mat<2> C(CMat,{nx,ny});

    // A= B+iC
    ComplexOperator<Mat<dim>,Mat<dim>,dim> A(B,C);

    GenEigsSolver<decltype(A)> eigs(A, 40, 65);

    // Initialize and compute
    eigs.init();
    int nconv = eigs.compute(SortRule::LargestMagn);

    // Retrieve results
    Eigen::VectorXcd evalues;
    if(eigs.info() == CompInfo::Successful)
        evalues = eigs.eigenvalues();

    auto error= evalues-evaluesOfA.template cast<std::complex<double>>();
    std::cout << error << std::endl;
    if (error.norm() > tol)
    {
        std::cout << "error= " << error.norm() << " is greater than tol " << tol << std::endl;
        return -1;

    }

    return 0;

}

