//
// Created by Nikhil Chandra Admal on 6/18/24.
//
#include <GbMaterialTensors.h>

namespace gbLAB {

    double GbMaterialTensors::tensorC(const int& k, const int& p, const int& l, const int& q)
    {
        return lambda*(k==p)*(l==q) + mu*(k==l)*(p==q) + mu*(k==q)*(p==l);
    }

    std::complex<double> GbMaterialTensors::tensorFhat(const int& k, const int& l, const int& i, const int& j, const Eigen::Vector3d& xi)
    {
        std::complex<double> output(0,0);
        double xi2Norm= sqrt(std::pow(xi(0),2) + std::pow(xi(2),2));
        double nu= lambda/(2*(lambda+mu));
        double nuFactor= 1-nu;
        // first term
        if (i==j)
        {
            if (k!=1 && l!=1)
                output= output + xi(k)*xi(l)*M_PI/(2*std::pow(xi2Norm,3));
            else if (k==1 && l==1)
                output= output + M_PI/(2*xi2Norm);
        }
        // second term
        std::list<int> ijkl{i,j,k,l};
        int count= std::count(ijkl.begin(),ijkl.end(),1);
        if ( count % 2 != 0)
            return output;
        else if(count == 0)
            output= output - xi(i)*xi(j)*xi(k)*xi(l)/(2*nuFactor) * 3 * M_PI / (8*std::pow(xi2Norm,5));
        else if(count == 4)
            output= output - 1.0/(2*nuFactor) * 3 * M_PI / (8*xi2Norm);
        else if(count == 2)
        {
            ijkl.remove(1);
            double temp=1.0;
            for(const auto elem : ijkl)
                temp= temp*xi(elem);
            output= output - temp/(2*nuFactor) * M_PI / (8*std::pow(xi2Norm,3));
        }

        return -output/(4*std::pow(M_PI,2)*mu);
    }


    std::complex<double> GbMaterialTensors::tensorGhat(const int& i, const int& k, const int& t, const int& r, const Eigen::VectorXd& xi)
    {
        std::complex<double> output(0,0);
        int l= 1;
        int s= 1;
        for (int u=0; u<3; ++u)
            for (int b=0; b<3; ++b)
                for (int c=0; c<3; ++c)
                    output= output + tensorC(i,l,u,r)*tensorC(b,c,t,s)*tensorFhat(c,k,u,b,xi) -
                            tensorC(i,l,u,s)*tensorC(b,c,t,r)*tensorFhat(c,k,u,b,xi) -
                            tensorC(i,k,u,r)*tensorC(b,c,t,s)*tensorFhat(c,l,u,b,xi) +
                            tensorC(i,k,u,s)*tensorC(b,c,t,r)*tensorFhat(c,l,u,b,xi);
        return output;
    }

    std::complex<double> GbMaterialTensors::tensorHhat(const int& t, const int& i, const Eigen::VectorXd &xi)
    {
        std::complex<double> output(0,0);
        for(int k=0; k<3; ++k)
            for(int r=0; r<3; ++r)
                output= output - 4*std::pow(M_PI,2)* (tensorGhat(i,k,t,r,xi) + tensorGhat(t,k,i,r,xi)) * xi(r)*xi(k);
        return output;
    }

    double GbMaterialTensors::lambda;
    double GbMaterialTensors::mu;
}

