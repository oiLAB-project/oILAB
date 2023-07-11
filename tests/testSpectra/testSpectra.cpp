#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Spectra/GenEigsSolver.h>
#include <iostream>
 
using namespace Spectra;

class MatProd
{
private:

    const Eigen::MatrixXd& M;

public:
    using Scalar = double;
    MatProd(const Eigen::MatrixXd& mat) : M(mat)
    {}
    Eigen::Index rows() const { return M.rows(); }
    Eigen::Index cols() const { return M.cols(); }
    void perform_op(const double* x_in, double* y_out) const
    {
            Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 1>> x(x_in, M.cols());
            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> y(y_out, M.rows());
            y.noalias() = M * x;
    }

};
 
int main()
{
    // A band matrix with 1 on the main diagonal, 2 on the below-main subdiagonal,
    // and 3 on the above-main subdiagonal
    const int n = 100;
    Eigen::Matrix<double,n,n> M;
    M.setZero();
    for(int i = 0; i < n; i++)
    {
        M(i, i) = 1.0;
        if(i > 0)
            M(i - 1, i) = 3.0;
        if(i < n - 1)
            M(i + 1, i) = 2.0;
    }

    // Construct matrix operation object using the wrapper class SparseGenMatProd
    //SparseGenMatProd<double> op(M);
    MatProd op(M);
    // Construct eigen solver object, requesting the largest three eigenvalues
    //GenEigsSolver<SparseGenMatProd<double>> eigs(op, 3, 6);
    GenEigsSolver<MatProd> eigs(op, 30, 65);

    // Initialize and compute
    eigs.init();
    int nconv = eigs.compute(SortRule::LargestMagn);
 
    // Retrieve results
    Eigen::VectorXcd evalues;
    if(eigs.info() == CompInfo::Successful)
        evalues = eigs.eigenvalues();
 
    std::cout << "Eigenvalues found:\n" << evalues << std::endl;
    std::cout << "\n Number of eigenvalues converged: " << nconv  << std::endl;

    return 0;
}

