//
// Created by Nikhil Chandra Admal on 6/23/23.
//

#ifndef OILAB_FFT_H
#define OILAB_FFT_H

#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <unsupported/Eigen/CXX11/Tensor>

class FFT : public Eigen::FFT<double>
{
public:
    using dcomplex= std::complex<double>;
    // 3 dimensions
    static void fft(const Eigen::Tensor<dcomplex,3>& in, Eigen::Tensor<dcomplex,3>& out)
    {
        out.setZero();
        Eigen::FFT<double> fft;
        for (int k = 0; k < in.dimension(0); k++)
        {
            for (int l = 0; l < in.dimension(1); l++)
            {
                Eigen::Tensor<dcomplex,1> tmpOutVector(in.dimension(2));
                tmpOutVector.setZero();
                FFT::fft(in.chip(k,0).chip(l,0),tmpOutVector);
                out.chip(k,0).chip(l,0)= tmpOutVector;
            }
        }

        for (int k = 0; k < in.dimension(1); k++)
        {
            for (int l = 0; l < in.dimension(2); l++)
            {
                Eigen::Tensor<dcomplex,1> tmpOutVector(in.dimension(0));
                tmpOutVector.setZero();
                FFT::fft(out.chip(k,1).chip(l,1),tmpOutVector);
                out.chip(k,1).chip(l,1)= tmpOutVector;
            }
        }

        for (int k = 0; k < in.dimension(0); k++)
        {
            for (int l = 0; l < in.dimension(2); l++)
            {
                Eigen::Tensor<dcomplex,1> tmpOutVector(in.dimension(0));
                tmpOutVector.setZero();
                FFT::fft(out.chip(k,0).chip(l,1),tmpOutVector);
                out.chip(k,0).chip(l,1)= tmpOutVector;
            }
        }
    }
    static void ifft(const Eigen::Tensor<dcomplex,3>& in, Eigen::Tensor<dcomplex,3>& out)
    {
        out.setZero();
        Eigen::FFT<double> fft;
        for (int k = 0; k < in.dimension(0); k++)
        {
            for (int l = 0; l < in.dimension(1); l++)
            {
                Eigen::Tensor<dcomplex,1> tmpOutVector(in.dimension(2));
                tmpOutVector.setZero();
                FFT::ifft(in.chip(k,0).chip(l,0),tmpOutVector);
                out.chip(k,0).chip(l,0)= tmpOutVector;
            }
        }

        for (int k = 0; k < in.dimension(1); k++)
        {
            for (int l = 0; l < in.dimension(2); l++)
            {
                Eigen::Tensor<dcomplex,1> tmpOutVector(in.dimension(0));
                tmpOutVector.setZero();
                FFT::ifft(out.chip(k,1).chip(l,1),tmpOutVector);
                out.chip(k,1).chip(l,1)= tmpOutVector;
            }
        }

        for (int k = 0; k < in.dimension(0); k++)
        {
            for (int l = 0; l < in.dimension(2); l++)
            {
                Eigen::Tensor<dcomplex,1> tmpOutVector(in.dimension(0));
                tmpOutVector.setZero();
                FFT::ifft(out.chip(k,0).chip(l,1),tmpOutVector);
                out.chip(k,0).chip(l,1)= tmpOutVector;
            }
        }
    }

    // 2 dimensions
    static void fft(const Eigen::Tensor<dcomplex,2>& in, Eigen::Tensor<dcomplex,2>& out)
    {
        out.setZero();
        Eigen::FFT<double> fft;

        for (int k = 0; k < in.dimension(0); k++) {
            Eigen::Tensor<dcomplex,1> tmpOutVector;
            tmpOutVector.setZero();
            FFT::fft(in.chip(k,0),tmpOutVector);
            out.chip(k,0)= tmpOutVector;
        }

        for (int k = 0; k < in.dimension(1); k++) {
            Eigen::Tensor<dcomplex,1> tmpOutVector;
            tmpOutVector.setZero();
            FFT::fft(out.chip(k,1),tmpOutVector);
            out.chip(k,1)= tmpOutVector;
        }
    }
    static void ifft(const Eigen::Tensor<dcomplex,2>& in, Eigen::Tensor<dcomplex,2>& out)
    {
        out.setZero();
        Eigen::FFT<double> fft;
        for (int k = 0; k < in.dimension(0); k++) {
            Eigen::Tensor<dcomplex,1> tmpOutVector(in.dimension(1));
            tmpOutVector.setZero();
            FFT::ifft(in.chip(k,0),tmpOutVector);
            out.chip(k,0)= tmpOutVector;
        }

        for (int k = 0; k < in.dimension(1); k++) {
            Eigen::Tensor<dcomplex,1> tmpOutVector(in.dimension(0));
            tmpOutVector.setZero();
            FFT::ifft(out.chip(k,1),tmpOutVector);
            out.chip(k,1)= tmpOutVector;
        }
    }

    // 1 dimension
    static void fft(const Eigen::Tensor<dcomplex,1>& in, Eigen::Tensor<dcomplex,1>& out)
    {
        out.setZero();
        Eigen::FFT<double> fft;
        Eigen::Map<const Eigen::VectorXcd> inMap(in.data(), in.dimension(0));

        Eigen::VectorXcd tmpOut(in.dimension(0));
        tmpOut.setZero();
        fft.fwd(tmpOut,inMap.eval());

        Eigen::TensorMap<Eigen::Tensor<dcomplex,1>> tmpOutMap(tmpOut.data(),tmpOut.size());
        out= tmpOutMap;
    }
    static void ifft(const Eigen::Tensor<dcomplex,1>& in, Eigen::Tensor<dcomplex,1>& out)
    {
        out.setZero();
        Eigen::FFT<double> fft;
        Eigen::Map<const Eigen::VectorXcd> inMap(in.data(), in.dimension(0));

        Eigen::VectorXcd tmpOut(in.dimension(0));
        tmpOut.setZero();
        fft.inv(tmpOut,inMap.eval());

        Eigen::TensorMap<Eigen::Tensor<dcomplex,1>> tmpOutMap(tmpOut.data(),tmpOut.size());
        out= tmpOutMap;
    }
};







// need to change in to a constant type


#endif //OILAB_FFT_H
