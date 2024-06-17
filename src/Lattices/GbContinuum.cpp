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
        bbhat(calculateb(domain, xuPairs, n)),
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
   std::pair<std::vector<PeriodicFunction<double,dim-1>>,
             std::vector<LatticeFunction<std::complex<double>,dim-1>>>
           GbContinuum<dim>::calculateb(const Eigen::Matrix<double,dim,dim-1>& domain,
                                         const std::deque<std::pair<VectorDimD,VectorDimD>>& xuPairs,
                                         const std::array<Eigen::Index,dim-1>& n)
   {
       Eigen::Matrix<double,dim,dim-1> basisVectors(domain.transpose().completeOrthogonalDecomposition().pseudoInverse());
       // lb0 - local b0
       std::vector<PeriodicFunction<double,dim-1>> lb0, lb;
       std::vector<LatticeFunction<std::complex<double>,dim-1>> lb0hat;
       for(int i=0; i<dim; ++i)
       {
           lb0.push_back(PeriodicFunction<double,dim-1>(n,domain));
           lb.push_back(PeriodicFunction<double,dim-1>(n,domain));
           //lb0hat.push_back(lb0[i].fft());
           lb0hat.push_back(LatticeFunction<std::complex<double>,dim-1>(n,basisVectors));
       }
       std::vector<LatticeFunction<std::complex<double>,dim-1>> lbhat(lb0hat);
       //Eigen::Matrix<double,dim,dim-1> basisVectors(lb0hat[0].basisVectors);



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
       Eigen::VectorXcd u2(numberOfCslPoints), lm(numberOfCslPoints);
       Eigen::MatrixXcd M(numberOfCslPoints,numberOfCslPoints);
       M.setZero();
       auto H22hatInverse= calculateH22hatInverse(domain,n);

       u2= u.col(2);
       for(int j=0; j<numberOfCslPoints; ++j)
       {
               for (int k=0; k<numberOfCslPoints; ++k)
               {
                   Eigen::Tensor<std::complex<double>, 0> sum;
                   Eigen::Tensor<std::complex<double>, dim - 1> temp;
                   temp= H22hatInverse.values * pihat[k].values * pihat[j].values.conjugate();
                   sum= temp.sum();
                   M(j,k)= sum(0);
               }
       }

       // Solve M lm = u2
       Eigen::LLT<Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic>> llt2;
       llt2.compute(M);
       lm= llt2.solve(u2);
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
       for(int k=0; k<numberOfCslPoints; ++k) {
           lbhat[2].values= lbhat[2].values + pihat[k].values * H22hatInverse.values * lm(k);
       }

       //lbhat[2].values= lbhat[2].values-lb0hat[2].values;

       for (int i=0; i<dim; ++i)
           lb[i].values = lbhat[i].ifft().values.real();

       //return std::make_pair(lb0,lb0hat);
       return std::make_pair(lb,lbhat);
   }

    template<int dim>
    GbContinuum<dim>::VectorDimD GbContinuum<dim>::displacement(const VectorDimD& t) const
    {
        //std::vector<LatticeFunction<std::complex<double>,dim-1>> b0hat;
        //for(int i=0; i<dim; ++i)
        //    b0hat.push_back(b0[i].fft());
        VectorDimD u;
        Eigen::Matrix<double,dim,dim-1> basisVectors(bhat[0].basisVectors);

        /*
        // function p
        VectorDimD normal(gbDomain.col(0).cross(gbDomain.col(1)));
        normal.normalize();
        DisplacementKernel f(normal,10);
        Shift<DisplacementKernel,double> p(t, (const Function<DisplacementKernel,double>&) f);
         */
        VectorDimD normal(gbDomain.col(0).cross(gbDomain.col(1)));
        normal.normalize();
        ShiftedDisplacementKernelFT phatFunction(t,normal);
        LatticeFunction<std::complex<double>,dim-1> phat(n,basisVectors,phatFunction);

        /*
         // debug
        LatticeFunction<std::complex<double>,dim-1> pphat(n,basisVectors);
        DisplacementKernel f(normal,13);
        Shift<DisplacementKernel,double> p(t, (const Function<DisplacementKernel,double>&) f);

        if (t.isZero()) {
            LatticeFunction<std::complex<double>, dim - 1> constantLatticeFunction(n, basisVectors);
            //phat.values.setConstant(1);
            pphat.values.setConstant(0.5);
        }
        else
            pphat.values= (p.fft<dim-1>(n, basisVectors)).values;

        std::cout << std::endl;
        std::cout << "t= " << t.transpose() << std::endl;
        std::cout << phat.values(1,0).real()/pphat.values(1,0).real() << std::endl;
        std::cout << phat.values(1,0).imag()/pphat.values(1,0).imag() << std::endl;
        std::cout << 1./(gbDomain.col(0).cross(gbDomain.col(1))).norm() << std::endl;
         */
        // debug

        /*
        // phat
        LatticeFunction<std::complex<double>,dim-1> phat(n,basisVectors);

        if (t.isZero()) {
            LatticeFunction<std::complex<double>, dim - 1> constantLatticeFunction(n, basisVectors);
            //phat.values.setConstant(1);
            phat.values.setConstant(0.5);
        }
        else
            phat.values= (p.fft<dim-1>(n, basisVectors)).values;

         */

        // u = f \star b
        for(int i=0; i<dim; ++i) {
            Eigen::Tensor<std::complex<double>, dim - 1> temp(bhat[i].values * phat.values.conjugate());
            Eigen::Tensor<std::complex<double>, 0> sum(temp.sum());
            u(i)= sum(0).real();
        }

        return u;

    }

    template<int dim>
    LatticeFunction<std::complex<double>,dim-1>
    GbContinuum<dim>::calculateH22hatInverse(const Eigen::Matrix<double, dim,dim-1>& domain,
                                      const std::array<Eigen::Index,dim-1>& n)
    {
        Eigen::Matrix<double,Eigen::Dynamic,dim-1> basisVectors(domain.transpose().completeOrthogonalDecomposition().pseudoInverse());

        H22InverseFT H22hatFunction(GbContinuum<dim>::lambda,GbContinuum<dim>::mu,basisVectors);
        return(LatticeFunction<std::complex<double>,dim-1>(n,basisVectors,H22hatFunction));
    }

    //template class GbContinuum<2>;
    template class GbContinuum<3>;
}
