//
// Created by Nikhil Chandra Admal on 5/27/24.
//
#include <GbContinuum.h>
#include <PeriodicFunction.cpp>
#include <Function.cpp>
#include <LatticeFunction.cpp>

namespace gbLAB {

    template<int dim>
    GbContinuum<dim>::GbContinuum(const Eigen::Matrix<double,dim,dim-1>& domain,
                                  const std::deque<std::pair<VectorDimD,VectorDimD>>& xuPairs,
                                  const std::array<Eigen::Index,dim-1>& n):
        gbDomain(domain),
        xuPairs(xuPairs),
        n(n),
        bbhat(calculateb(domain,xuPairs,n)),
        b(bbhat.first),
        bhat(bbhat.second)
   {
       //check

       /*
       for(const auto& xu : xuPairs) {
           std::cout << "u = " << xu.second.transpose() << std::endl;
           std::cout << "u solved = " << displacement(xu.first).transpose() << std::endl;
       }
       std::cout << std::endl;
        */
   }

   template<int dim>
   typename GbContinuum<dim>::FunctionFFTPair
   GbContinuum<dim>::calculateb(const Eigen::Matrix<double,dim,dim-1>& domain,
                                const std::deque<std::pair<VectorDimD,VectorDimD>>& xuPairs,
                                const std::array<Eigen::Index,dim-1>& n)
   {
       if (HhatInvComponents.size() ==0 ) HhatInvComponents= getHhatInvComponents(domain,n);
       Eigen::Matrix<double,dim,dim-1> basisVectors(domain.transpose().completeOrthogonalDecomposition().pseudoInverse());
       // lb0 - read as "local b0"
       std::vector<PeriodicFunction<double,dim-1>> lb0, lb;
       std::vector<LatticeFunction<std::complex<double>,dim-1>> lb0hat;
       for(int i=0; i<dim; ++i)
       {
           lb0.push_back(PeriodicFunction<double,dim-1>(n,domain));
           lb.push_back(PeriodicFunction<double,dim-1>(n,domain));
           lb0hat.push_back(LatticeFunction<std::complex<double>,dim-1>(n,basisVectors));
       }
       std::vector<LatticeFunction<std::complex<double>,dim-1>> lbhat(lb0hat);


       // calculate GB normal
       VectorDimD normal(domain.col(0).cross(domain.col(1)));
       normal.normalize();

       // calculare pi and pihat
       std::vector<PeriodicFunction<double,dim-1>> pi;
       std::vector<LatticeFunction<std::complex<double>,dim-1>> pihat;
       DisplacementKernel f(normal,10);

       for(const auto& xu : xuPairs) {
           VectorDimD t = xu.first;
           ShiftedDisplacementKernelFT pihatFunction(t,normal);
           Shift<DisplacementKernel,double> piFunction(t, (const Function<DisplacementKernel,double>&) f);
           pi.push_back(PeriodicFunction<double,dim-1>(n,domain,piFunction));
           pihat.push_back(LatticeFunction<std::complex<double>,dim-1> (n,basisVectors,pihatFunction));
           auto temp= *(pihat.end()-1);
       }

       // matrices P and u
       const int numberOfCslPoints= xuPairs.size();
       Eigen::Matrix<std::complex<double>,Eigen::Dynamic,dim> u(numberOfCslPoints,dim);
       Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic> P(numberOfCslPoints,numberOfCslPoints);
       for(int i=0; i<numberOfCslPoints; ++i)
       {
           u.row(i) = xuPairs[i].second;
           for (int j=0; j<numberOfCslPoints; ++j)
           {
               Eigen::Tensor<std::complex<double>,dim-1> temp(pihat[j].values * pihat[i].values.conjugate());
               Eigen::Tensor<std::complex<double>,0> sum(temp.sum());
               P(i,j)= sum(0);
           }
       }

       // Solve P \alpha = u
       Eigen::Matrix<std::complex<double>,Eigen::Dynamic,dim> alpha(numberOfCslPoints,dim);
       Eigen::LLT<Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic>> llt;
       llt.compute(P);
       alpha= llt.solve(u);

       // calculate lb0 and lb0hat
       for(int i=0; i<dim; ++i)
       {
           lb0hat[i].values.setZero();
           for(int j=0; j<numberOfCslPoints; ++j)
               lb0hat[i].values= lb0hat[i].values + pihat[j].values * alpha(j,i);
       }
       for(int i=0; i<dim; ++i)
           for(int j=0; j<numberOfCslPoints; ++j)
               lb0[i].values= lb0[i].values + pi[j].values * alpha(j,i).real();


       // calculate lagrange multipliers, lbhat, and lb
       Eigen::VectorXcd uFlattened;
       uFlattened= u.reshaped();
       Eigen::VectorXcd lm(dim*numberOfCslPoints);
       Eigen::MatrixXcd M(dim*numberOfCslPoints,dim*numberOfCslPoints);
       M.setZero();

       for (int i=0; i<dim; ++i)
       {
           for (int j=0; j<numberOfCslPoints; ++j)
           {
               int row= i*numberOfCslPoints + j;
               for (int l=0; l<dim; ++l)
               {
                   for (int k = 0; k < numberOfCslPoints; ++k)
                   {
                       int col= l*numberOfCslPoints + k;
                       Eigen::Tensor<std::complex<double>, 0> sum;
                       Eigen::Tensor<std::complex<double>, dim - 1> temp(n);
                       LatticeFunction<std::complex<double>,dim-1> HhatInv_il(n,basisVectors);
                       if (i==0 && l==0) temp = HhatInvComponents[0].values * pihat[k].values * pihat[j].values.conjugate();
                       if (i==1 && l==1) temp = HhatInvComponents[1].values * pihat[k].values * pihat[j].values.conjugate();
                       if (i==2 && l==2) temp = HhatInvComponents[2].values * pihat[k].values * pihat[j].values.conjugate();
                       if ((i==1 && l==2) || (i==2 && l==1)) temp = HhatInvComponents[3].values * pihat[k].values * pihat[j].values.conjugate();
                       if ((i==0 && l==2) || (i==2 && l==0)) temp = HhatInvComponents[4].values * pihat[k].values * pihat[j].values.conjugate();
                       if ((i==0 && l==1) || (i==1 && l==0)) temp = HhatInvComponents[5].values * pihat[k].values * pihat[j].values.conjugate();
                       sum = temp.sum();
                       M(row, col) = sum(0);
                   }
               }
           }
       }

       // Solve M lm = u2
       Eigen::LLT<Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic>> llt2;
       llt2.compute(M);
       lm= llt2.solve(uFlattened);
       /*
       std::cout << "u = " << u.col(2).transpose() << std::endl;
       std::cout << "M = " << M << std::endl;
       std::cout << "P = " << P << std::endl;
       std::cout << "lm = " << lm.transpose() << std::endl;
       std::cout << "alpha = " << alpha << std::endl;
        */

       lbhat[0].values.setZero();
       lbhat[1].values.setZero();
       lbhat[2].values.setZero();

       auto lmMatrix= lm.reshaped(numberOfCslPoints,dim);
       // temp = \sum_k p_k * lm_k
       std::vector<LatticeFunction<std::complex<double>,dim-1>> temp;

       for (int i=0; i<dim; ++i) {
           temp.push_back(LatticeFunction<std::complex<double>, dim - 1>(n, basisVectors));
           for (int k = 0; k < numberOfCslPoints; ++k)
               temp[i].values = temp[i].values + pihat[k].values * lmMatrix(k, i);
       }

       lbhat[0].values= HhatInvComponents[0].values * temp[0].values +
                        HhatInvComponents[5].values * temp[1].values +
                        HhatInvComponents[4].values * temp[2].values;

       lbhat[1].values= HhatInvComponents[5].values * temp[0].values +
                        HhatInvComponents[1].values * temp[1].values +
                        HhatInvComponents[3].values * temp[2].values;

       lbhat[2].values= HhatInvComponents[4].values * temp[0].values +
                        HhatInvComponents[3].values * temp[1].values +
                        HhatInvComponents[2].values * temp[2].values;

       for (int i=0; i<dim; ++i)
           lb[i].values = lbhat[i].ifft().values.real();

       return std::make_pair(lb,lbhat);
   }

    template<int dim>
    GbContinuum<dim>::VectorDimD GbContinuum<dim>::displacement(const VectorDimD& t) const
    {
        VectorDimD u;
        Eigen::Matrix<double,dim,dim-1> basisVectors(bhat[0].basisVectors);

        VectorDimD normal(gbDomain.col(0).cross(gbDomain.col(1)));
        normal.normalize();
        ShiftedDisplacementKernelFT phatFunction(t,normal);
        LatticeFunction<std::complex<double>,dim-1> phat(n,basisVectors,phatFunction);

        // u = f \star b
        for(int i=0; i<dim; ++i) {
            Eigen::Tensor<std::complex<double>, dim - 1> temp(bhat[i].values * phat.values.conjugate());
            Eigen::Tensor<std::complex<double>, 0> sum(temp.sum());
            u(i)= sum(0).real();
        }

        return u;

    }

    template<int dim>
    std::vector<LatticeFunction<std::complex<double>,dim-1>>
    GbContinuum<dim>::getHhatInvComponents(const Eigen::Matrix<double, dim,dim-1>& domain,
                                           const std::array<Eigen::Index,dim-1>& n)
    {
        Eigen::Matrix<double,Eigen::Dynamic,dim-1> basisVectors(domain.transpose().completeOrthogonalDecomposition().pseudoInverse());
        std::vector<LatticeFunction<std::complex<double>,dim-1>> output;

        output.push_back(LatticeFunction<std::complex<double>,dim-1>(n,basisVectors,HhatInvFunction(0,0,domain)));
        output.push_back(LatticeFunction<std::complex<double>,dim-1>(n,basisVectors,HhatInvFunction(1,1,domain)));
        output.push_back(LatticeFunction<std::complex<double>,dim-1>(n,basisVectors,HhatInvFunction(2,2,domain)));
        output.push_back(LatticeFunction<std::complex<double>,dim-1>(n,basisVectors,HhatInvFunction(1,2,domain)));
        output.push_back(LatticeFunction<std::complex<double>,dim-1>(n,basisVectors,HhatInvFunction(0,2,domain)));
        output.push_back(LatticeFunction<std::complex<double>,dim-1>(n,basisVectors,HhatInvFunction(0,1,domain)));
        return output;
    }


    /*----------------------------*/
    DisplacementKernel::DisplacementKernel(Eigen::Vector<double, Eigen::Dynamic> _normal,
                                           double domainSize) :
        normal(_normal), Function<DisplacementKernel, double>(domainSize)
    {}
    double DisplacementKernel::operator()(const Eigen::Vector<double, Eigen::Dynamic> &x) const
    {
        return -x.dot(normal) / (4 * M_PI * (std::pow(x.norm(), 3)));
    }
    /*----------------------------*/
    ShiftedDisplacementKernelFT::ShiftedDisplacementKernelFT(const Eigen::VectorXd& x,
                                                             const Eigen::VectorXd& normal):
            x(x),normal(normal.normalized())
    { }

    std::complex<double> ShiftedDisplacementKernelFT::operator()(const Eigen::VectorXd& xi) const
    {
        double xNormalComponent=  x.dot(normal);
        auto xProjected= (x-xNormalComponent*normal).eval();
        double xiBar= xi.norm();
        std::complex<double> output= -0.5*std::exp(-2*M_PI* std::complex<double>(0,1) *
                                                   xProjected.dot(xi)
        ) * std::exp(-2*M_PI*xiBar*abs(xNormalComponent));
        if (abs(xNormalComponent) < DBL_EPSILON)
            return output;
        else
            return output*xNormalComponent/abs(xNormalComponent);
    }
    /*----------------------------*/
    HhatInvFunction::HhatInvFunction(const int& t, const int& i, const Eigen::MatrixXd& domain) :
            t(t),i(i), e1(domain.col(0).normalized()), e3(domain.col(1).normalized())
    { }
    std::complex<double> HhatInvFunction::operator()(const Eigen::VectorXd& xi) const
    {
        if (xi.isZero())
            return std::complex<double>(0,0);
        std::complex<double> output(0,0);
        double xi1= xi.dot(e1);
        double xi3= xi.dot(e3);
        double xi2= (xi-xi1*e1-xi3*e3).norm();

        output= tensorHhat(t,i,(Eigen::VectorXd(3) << xi1,xi2,xi3).finished());
        output= -std::pow(2*M_PI,2) * output;
        return output;
    }
    /*----------------------------*/

    // Explicit instantiation
    //template class GbContinuum<2>;
    template class GbContinuum<3>;
}
