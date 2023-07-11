#include<iostream>
#include <Eigen/Core>
#include<LatticeModule.h>
#include<FFT.h>

// A test that computes the Laplacian using FFT
using namespace gbLAB;
int main ()
{
    using dcomplex=std::complex<double>;
    double tol=1e-8;
    Eigen::Matrix2d A;
    //A << 1.0 , 0.0,
    //        0.0, 1.0;
    A << 11.2 , 3.5,
         2.0, 15.0;

    Lattice<2> L(A);
    LatticeVector<2> l1(L);
    LatticeVector<2> l2(L);
    l1 << 1,0; l2 << 0,1;

    ReciprocalLatticeVector<2> rl1(L);
    ReciprocalLatticeVector<2> rl2(L);
    rl1 << 1,0; rl2 << 0,1;

    int Nx = 32, Ny= 32;
    //Eigen::MatrixXd f(Nx,Ny);
    Eigen::Tensor<dcomplex,2> f(Nx,Ny);
    f.setZero();

    for(int i=0; i<Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            Eigen::Vector2d l;
            l << (double) i * l1.cartesian() / Nx + (double) j * l2.cartesian() / Ny;
            // f = sin(2 pi x.r1)  cos(2 pi x.r2)
            f(i, j) = cos(2*M_PI*l.dot(rl1.cartesian())) *  cos(2*M_PI*l.dot(rl2.cartesian()));
        }
    }


    // compute Laplacian
    Eigen::MatrixXcd fhat;
    fhat.setZero(Nx, Ny);
    FFT transform;
    transform.fft(f,fhat);

    Eigen::MatrixXcd d2fhat;
    d2fhat.setZero(Nx, Ny);
    for(int i=0; i<Nx; ++i)
    {
        for (int j = 0; j < Ny; ++j) {
            ReciprocalLatticeVector<2> r(L);
            //r << (i == Nx / 2 ? 0 : -Nx * (i/(Nx / 2)) + i), (j == Ny / 2 ? 0 : -Ny * (j / (Ny / 2)) + j);
            r << (i <= Nx / 2 ? i : i-Nx ), (j <= Ny / 2 ? j : j-Ny );
            d2fhat(i, j) = -fhat(i, j) * (std::complex<double>) (4.0 * std::pow(M_PI, 2) * r.cartesian().squaredNorm());
        }
    }

    // Laplacian Lf
    Eigen::MatrixXcd Lf;
    Lf.setZero(Nx, Ny);
    transform.ifft(d2fhat,Lf);

    // Analytical Laplacian aLf = -4 pi^2 (|r1|^2 + |r2|^2) sin(2 pi x.r1)  cos(2 pi x.r2) - 2 r1.r2 cos(2 pi x.r1) sin(2 pi x.r2)
    Eigen::MatrixXcd aLf(Nx,Ny);
    Eigen::MatrixXcd adxf(Nx,Ny);
    aLf.setZero();
    adxf.setZero();
    for(int i=0; i<Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            Eigen::Vector2d l;
            l << (double) i * l1.cartesian() / Nx + (double) j * l2.cartesian() / Ny;
            aLf(i, j) = -rl1.cartesian().squaredNorm()*cos(2*M_PI*l.dot(rl1.cartesian())) *  cos(2*M_PI*l.dot(rl2.cartesian()))
                                 -rl2.cartesian().squaredNorm()*cos(2*M_PI*l.dot(rl1.cartesian())) *  cos(2*M_PI*l.dot(rl2.cartesian()))
                                 +2*rl1.cartesian().dot(rl2.cartesian())*sin(2*M_PI*l.dot(rl1.cartesian())) *  sin(2*M_PI*l.dot(rl2.cartesian()));
            adxf(i, j) = -rl1.cartesian().x()*sin(2*M_PI*l.dot(rl1.cartesian())) *  cos(2*M_PI*l.dot(rl2.cartesian()))
                                 -rl2.cartesian().x()*cos(2*M_PI*l.dot(rl1.cartesian())) *  sin(2*M_PI*l.dot(rl2.cartesian()));
        }
    }
    aLf=4*std::pow(M_PI,2)*aLf;
    adxf=2*M_PI*adxf;

    std::cout << "f= " << std::endl;
    std::cout << f.real() << std::endl;
    std::cout << "Lf = " << std::endl;
    std::cout << Lf.real() << std::endl;

    std::cout << "Analytical Lf = " << std::endl;
    std::cout << aLf.real() << std::endl;

    double error= (Lf.real()-aLf.real()).norm()/(Nx*Ny);
    std::cout << "error in Laplacian= " << std::endl;
    std::cout << error << std::endl;
    if (error >= tol)
        throw std::runtime_error("Test FFT failed");


    Eigen::MatrixXcd dxfhat;
    dxfhat.setZero(Nx, Ny);
    for(int i=0; i<Nx; ++i)
    {
        for (int j = 0; j < Ny; ++j) {
            ReciprocalLatticeVector<2> r(L);
            r << (i == Nx / 2 ? 0 : -Nx * (i/(Nx / 2)) + i), (j == Ny / 2 ? 0 : -Ny * (j / (Ny / 2)) + j);
            //r << (i <= Nx / 2 ? i : i-Nx ), (j <= Ny / 2 ? j : j-Ny );
            dxfhat(i, j) = fhat(i, j) * (std::complex<double>) (2.0 * M_PI * 1i * r.cartesian().x());
        }
    }
    Eigen::MatrixXcd dxf;
    dxf.setZero(Nx, Ny);
    transform.ifft(dxfhat,dxf);

    std::cout << adxf << std::endl;

    std::cout << "error in dxf = " << std::endl;
    error= (dxf.real()-adxf.real()).norm()/(Nx*Ny);
    std::cout <<  error << std::endl;
    if (error >= tol)
        throw std::runtime_error("Test FFT failed");


    return 0;
}
