#include<iostream>
#include<LatticeModule.h>
#include<fstream>
#include<FFT.h>

#define MY_ERROR(message)                                                \
  {                                                                      \
    std::cout << "* Error : \"" << message << "\" : " << __LINE__ << ":" \
              << __FILE__ << std::endl;                                  \
    exit(1);                                                             \
  }

// A test that computes the Laplacian using FFT
using namespace gbLAB;
int main ()
{
    using dcomplex= std::complex<double>;
    std::ifstream file3D("m3D.txt");
    if(!file3D) MY_ERROR(std::string("ERROR: m3D.txt could not be opened for reading!"));
    std::ifstream file2D("m2D.txt");
    if(!file2D) MY_ERROR(std::string("ERROR: m2D.txt could not be opened for reading!"));

    std::ifstream fileFFT3D("fft_m3D.txt");
    if(!fileFFT3D) MY_ERROR(std::string("ERROR: fft_m3D.txt could not be opened for reading!"));
    std::ifstream fileFFT2D("fft_m2D.txt");
    if(!fileFFT2D) MY_ERROR(std::string("ERROR: fft_m2D.txt could not be opened for reading!"));

    double tol= 1e-6;

    Eigen::Tensor<dcomplex,3> m3D(10,15,20);
    Eigen::Tensor<dcomplex,3> n3DRef(10,15,20);

    Eigen::Tensor<dcomplex,2> m2D(15,20);
    Eigen::Tensor<dcomplex,2> n2DRef(15,20);

    for (int k=0; k<20; k++)
        for (int j=0; j<15; j++)
        {
            file2D >> m2D(j,k);
            double tmp1, tmp2;
            fileFFT2D >> tmp1 >> tmp2;
            n2DRef(j, k) =dcomplex (tmp1, tmp2);
            for (int i = 0; i < 10; i++)
            {
                file3D >> m3D(i, j, k);
                fileFFT3D >> tmp1 >> tmp2;
                n3DRef(i, j, k) =dcomplex(tmp1,tmp2);
            }
        }

    Eigen::Tensor<dcomplex,3> n3D(10,15,20);
    Eigen::Tensor<dcomplex,2> n2D(15,20);

    FFT::fft(m3D,n3D);
    FFT::fft(m2D,n2D);

    // check error
    Eigen::Tensor<double,0> error3D= (n3D-n3DRef).abs().maximum(Eigen::array<int, 3>({0, 1,2}) );
    std::cout << "error in forward 3D fft = " << error3D(0) << std::endl;
    if (error3D(0) > tol)
        return -1;

    Eigen::Tensor<double,0> error2D= (n2D-n2DRef).abs().maximum(Eigen::array<int, 2>({0, 1}) );
    std::cout << "error in forward 2D fft = " << error2D(0) << std::endl;
    if (error2D(0) > tol)
        return -1;


    Eigen::Tensor<dcomplex,3> m3DRecovered(10,15,20);
    FFT::ifft(n3D,m3DRecovered);
    error3D= (m3D.cast<dcomplex>()-m3DRecovered).abs().maximum(Eigen::array<int, 3>({0, 1,2}) );
    std::cout << "error in inverse 3D fft = " << error3D(0) << std::endl;
    if (error3D(0) > tol)
        return -1;

    Eigen::Tensor<dcomplex,2> m2DRecovered(15,20);
    FFT::ifft(n2D,m2DRecovered);
    //error2D= (m2D.cast<dcomplex>()-m2DRecovered).norm();
    error2D= (m2D.cast<dcomplex>()-m2DRecovered).abs().maximum(Eigen::array<int, 2>({0, 1}) );
    std::cout << "error in inverse 2D fft = " << error2D(0) << std::endl;
    if (error2D(0) > tol)
        return -1;
    return 0;
}
