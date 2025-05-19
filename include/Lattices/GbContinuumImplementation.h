//
// Created by Nikhil Chandra Admal on 5/27/24.
//

#ifndef OILAB_GBCONTINUUMIMPLEMENTATION_H
#define OILAB_GBCONTINUUMIMPLEMENTATION_H

#include<Diff.h>

namespace gbLAB {

    template<int dim>
    GbContinuum<dim>::GbContinuum(const Eigen::Matrix<double,dim,dim-1>& domain,
                                  const std::map<OrderedTuplet<dim+1>,VectorDimD>& xuPairs,
                                  const std::array<Eigen::Index,dim-1>& n,
                                  const std::map<OrderedTuplet<dim+1>,VectorDimD>& atoms,
                                  const bool& verbosity):
        gbDomain(domain),
        xuPairs(xuPairs),
        n(n),
        bbhat(calculateb(domain,xuPairs,n,atoms)),
        b(bbhat.first),
        bhat(bbhat.second),
        alpha(get_alpha(b))
   {
       uAverage.setZero();
       for (const auto &[x, u]: xuPairs) {
           uAverage = uAverage + u;
       }
       if (xuPairs.size() != 0)
           uAverage = uAverage / xuPairs.size();

       if (verbosity) {
           std::cout << "------------------------------------------------------------------------------" << std::endl;
           std::cout << "Constraints: " << std::endl;
       }
       for (const auto &[x, u]: xuPairs) {
           if (verbosity)
               std::cout << "x = " << atoms.at(x).transpose() << "; displacement = " << u.transpose() << std::endl;
           if ((u - displacement(x) - uAverage).norm() > FLT_EPSILON)
               throw std::runtime_error("GBContinuum construction failed - unable to impose constraints.");
       }
       if (verbosity)
           std::cout << std::endl;

   }

    template<int dim>
    std::vector<PeriodicFunction<double,dim-1>>
    GbContinuum<dim>::get_alpha(const std::vector<PeriodicFunction<double, dim - 1>>& b)
    {
        Eigen::Matrix<double,dim,dim-1> orthogonalVectors;
        Eigen::Matrix<double,dim,dim-1> unitCell= b[0].unitCell;
        Eigen::array<Eigen::Index,dim-1> n= b[0].values.dimensions();

        // Form a rotation matrix to transform to a 2D coordinate system containing the grain boundary
        for(int i=0; i<dim-1; ++i)
        {
            orthogonalVectors.col(i)= unitCell.col(i);
            for(int j=0; j<i; ++j)
                orthogonalVectors.col(i)-= unitCell.col(i).dot(orthogonalVectors.col(j))* orthogonalVectors.col(j);
            orthogonalVectors.col(i)= orthogonalVectors.col(i).normalized();
        }
        Rotation<dim> rotation(orthogonalVectors);
        assert((rotation*unitCell).row(2).norm()<FLT_EPSILON);

        // transform the slip vectors b to the local coordinate system
        std::vector<PeriodicFunction<double, dim - 1>>  lb;
        for(int i= 0; i<dim; ++i) {
            PeriodicFunction<double, dim - 1> temp(n,unitCell);
            for (int j = 0; j < dim; ++j)
                temp.values = temp.values + rotation(i, j) * b[j].values;
            lb.push_back(temp);
        }

        // Define the gradient operator
        Eigen::Matrix<double,dim-1,dim-1> A= (rotation*unitCell).block(0,0,dim-1,dim-1);
        // dimension-dependent part
        Diff<dim-1> dx({1,0},A,n);
        Diff<dim-1> dy({0,1},A,n);

        std::vector<PeriodicFunction<double,dim-1>> alpha;
        for(int i=0; i<dim; ++i)
        {
            for(int j=0; j<dim-1; ++j)
            {
                PeriodicFunction<double,dim-1> alphaij(n, unitCell);;
                for(int k=0; k<dim-1; ++k)
                {
                    if (k==j) continue;
                    PeriodicFunction<double,dim-1> dlbi_dxk(n, unitCell);;
                    if (k==0)
                        dx.perform_op(lb[i].values.data(), dlbi_dxk.values.data());
                    else if (k==1)
                        dy.perform_op(lb[i].values.data(), dlbi_dxk.values.data());
                    if (j==0 && k==1)
                        alphaij.values= alphaij.values + dlbi_dxk.values;
                    else if (j==1 && k==0)
                        alphaij.values= alphaij.values - dlbi_dxk.values;
                }
                alpha.push_back(alphaij);
            }
        }

        return alpha;
    }

    template<int dim>
    PeriodicFunction<double,dim-1>
    GbContinuum<dim>::get_pi(const Eigen::Matrix<double,dim,dim-1>& domain,
                             const std::array<Eigen::Index,dim-1>& n,
                             const VectorDimD& point)
    {
        Eigen::Matrix<double,dim,dim-1> basisVectors(domain.transpose().completeOrthogonalDecomposition().pseudoInverse());
        VectorDimD normal(domain.col(0).cross(domain.col(1)));
        normal.normalize();

        //DisplacementKernel f(normal,10);
        // change 1e6 entry
        DisplacementKernel f(normal);
        Shift<DisplacementKernel,double> piFunction(point, (const Function<DisplacementKernel,double>&) f);
        return PeriodicFunction<double,dim-1>(n,domain,piFunction);
    }


    template<int dim>
    LatticeFunction<std::complex<double>,dim-1>
    GbContinuum<dim>::get_pihat(const Eigen::Matrix<double,dim,dim-1>& domain,
                                const std::array<Eigen::Index,dim-1>& n,
                                const VectorDimD& point)
    {
        Eigen::Matrix<double,dim,dim-1> basisVectors(domain.transpose().completeOrthogonalDecomposition().pseudoInverse());
        VectorDimD normal(domain.col(0).cross(domain.col(1)));
        normal.normalize();

        std::vector<LatticeFunction<std::complex<double>,dim-1>> pihat;
        //DisplacementKernel f(normal,10);

        ShiftedDisplacementKernelFT pihatFunction(point,normal);
        return LatticeFunction<std::complex<double>,dim-1> (n,basisVectors,pihatFunction);
    }

   template<int dim>
   typename GbContinuum<dim>::FunctionFFTPair
   GbContinuum<dim>::calculateb(const Eigen::Matrix<double,dim,dim-1>& domain,
                                const std::map<OrderedTuplet<dim+1>,VectorDimD>& xuPairs,
                                const std::array<Eigen::Index,dim-1>& n,
                                const std::map<OrderedTuplet<dim+1>,VectorDimD>& atoms)
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


       if(piPeriodicFunctions.empty()) {
           for (const auto& [key, value]: atoms) {
               // the cross product of the domain vectors has to be parallel to nA
               Eigen::Matrix<double,dim,dim-1> basisVectors(domain.transpose().completeOrthogonalDecomposition().pseudoInverse());
               VectorDimD normal(domain.col(0).cross(domain.col(1)));
               normal.normalize();

               VectorDimD perturbedValue= value;
               if(abs(value.dot(normal)) < FLT_EPSILON) // belongs to the GB
               {
                   if (key(dim)==1)
                       perturbedValue= value - (value.dot(normal) + FLT_EPSILON) * normal;
                   else if (key(dim)==2)
                       perturbedValue= value - (value.dot(normal) - FLT_EPSILON) * normal;
               }

               else if (key(dim)==-1 || key(dim)==-2) // belongs to (lattice 1 & grain 2) || (lattice 2 & grain 1)
               {
                   perturbedValue= value - 2 * value.dot(normal) * normal;
               }
               piPeriodicFunctions.insert({key, get_pi(domain, n, perturbedValue)});
               pihatLatticeFunctions.insert({key, get_pihat(domain, n, perturbedValue)});

               /*
               if(abs(value.dot(normal)) < FLT_EPSILON && key(dim)==1) // belongs to lattice 1 and on the GB
               {
                   VectorDimD perturbedValue= value - (value.dot(normal) + FLT_EPSILON) * normal;
                   piPeriodicFunctions.insert({key, get_pi(domain, n, perturbedValue)});
                   pihatLatticeFunctions.insert({key, get_pihat(domain, n, perturbedValue)});
               }
               else if(abs(value.dot(normal)) < FLT_EPSILON && key(dim)==2) // belongs to lattice 2 and on the GB
               {
                   VectorDimD perturbedValue= value - (value.dot(normal) - FLT_EPSILON) * normal;
                   piPeriodicFunctions.insert({key, get_pi(domain, n, perturbedValue)});
                   pihatLatticeFunctions.insert({key, get_pihat(domain, n, perturbedValue)});
               }
               else {
                   piPeriodicFunctions.insert({key, get_pi(domain, n, value)});
                   pihatLatticeFunctions.insert({key, get_pihat(domain, n, value)});
               }
                */
           }
       }

       // matrices P and u
       const int numberOfCslPoints= xuPairs.size();
       Eigen::Matrix<std::complex<double>,Eigen::Dynamic,dim> uMatrix(numberOfCslPoints,dim);
       Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic> P(numberOfCslPoints,numberOfCslPoints);


       VectorDimD uAverage;
       uAverage.setZero();
       for(const auto& [xi,ui] : xuPairs)
       {
           uAverage= uAverage + ui;
       }
       if (xuPairs.size() != 0)
           uAverage= uAverage/xuPairs.size();

       int row= -1;
       for(const auto& [xi,ui] : xuPairs)
       {
           row++;
           //uMatrix.row(row)= ui;
           uMatrix.row(row)= ui-uAverage;
           int col= -1;
           for(const auto& [xj,uj] : xuPairs)
           {
               col++;
               P(row, col) = pihatLatticeFunctions.at(xj).dot(pihatLatticeFunctions.at(xi));
           }

       }

       // Solve P \alpha = u
       Eigen::Matrix<std::complex<double>,Eigen::Dynamic,dim> alpha(numberOfCslPoints,dim);
       Eigen::LLT<Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic>> llt;
       llt.compute(P);
       alpha= llt.solve(uMatrix);

       // calculate lb0 and lb0hat

       for(int i=0; i<dim; ++i)
       {
           lb0hat[i].values.setZero();
           int j= -1;
           for(const auto& [xj,uj] : xuPairs) {
               j++;
               lb0hat[i].values = lb0hat[i].values + pihatLatticeFunctions.at(xj).values * alpha(j,i);
           }
       }
       for(int i=0; i<dim; ++i) {
           lb0[i].values.setZero();
           int j = -1;
           for (const auto &[xj, uj]: xuPairs) {
               j++;
               lb0[i].values = lb0[i].values + piPeriodicFunctions.at(xj).values * alpha(j,i).real();
           }
       }

       // calculate lagrange multipliers, lbhat, and lb
       Eigen::VectorXcd uFlattened;
       uFlattened= uMatrix.reshaped();
       Eigen::VectorXcd lm(dim*numberOfCslPoints);
       Eigen::MatrixXcd M(dim*numberOfCslPoints,dim*numberOfCslPoints);
       M.setZero();

       for (int i=0; i<dim; ++i)
       {
           int j= -1;
           for (const auto& [xj,uj] : xuPairs)
           {
               ++j;
               int row= i*numberOfCslPoints + j;
               for (int l=0; l<dim; ++l)
               {
                   int k= -1;
                   for(const auto& [xk,uk] : xuPairs)
                   {
                       k++;
                       int col= l*numberOfCslPoints + k;

                       std::complex<double> temp;

                       if (i==0 && l==0) temp = (HhatInvComponents[0]*pihatLatticeFunctions.at(xk)).dot(pihatLatticeFunctions.at(xj));
                       if (i==1 && l==1) temp = (HhatInvComponents[1]*pihatLatticeFunctions.at(xk)).dot(pihatLatticeFunctions.at(xj));
                       if (i==2 && l==2) temp = (HhatInvComponents[2]*pihatLatticeFunctions.at(xk)).dot(pihatLatticeFunctions.at(xj));
                       if ((i==1 && l==2) || (i==2 && l==1)) temp = (HhatInvComponents[3]*pihatLatticeFunctions.at(xk)).dot(pihatLatticeFunctions.at(xj));
                       if ((i==0 && l==2) || (i==2 && l==0)) temp = (HhatInvComponents[4]*pihatLatticeFunctions.at(xk)).dot(pihatLatticeFunctions.at(xj));
                       if ((i==0 && l==1) || (i==1 && l==0)) temp = (HhatInvComponents[5]*pihatLatticeFunctions.at(xk)).dot(pihatLatticeFunctions.at(xj));
                       M(row, col) = temp;
                   }
               }
           }
       }

       // Solve M lm = u2
       Eigen::LLT<Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic>> llt2;
       llt2.compute(M);
       lm= llt2.solve(uFlattened);

       lbhat[0].values.setZero();
       lbhat[1].values.setZero();
       lbhat[2].values.setZero();

       auto lmMatrix= lm.reshaped(numberOfCslPoints,dim);
       // temp = \sum_k p_k * lm_k
       std::vector<LatticeFunction<std::complex<double>,dim-1>> temp;

       for (int i=0; i<dim; ++i) {
           temp.push_back(LatticeFunction<std::complex<double>, dim - 1>(n, basisVectors));
           int k= -1;
           for (const auto& [xk,uk] : xuPairs) {
               k++;
               temp[i].values = temp[i].values + pihatLatticeFunctions.at(xk).values * lmMatrix(k, i);
           }
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

       //return std::make_pair(lb,lbhat);
       return std::make_pair(lb,lbhat);
   }

    template<int dim>
    typename GbContinuum<dim>::VectorDimD GbContinuum<dim>::displacement(const OrderedTuplet<dim+1>& t) const
    {
        VectorDimD u;

        // u = f \star b
        for(int i=0; i<dim; ++i)
            u(i)= (bhat[i].dot(pihatLatticeFunctions.at(t))).real();

        return u;

    }


    template<int dim>
    typename GbContinuum<dim>::VectorDimD GbContinuum<dim>::displacement(const VectorDimD& t) const
    {
        VectorDimD u;

        Eigen::Matrix<double,dim,dim-1> basisVectors(gbDomain.transpose().completeOrthogonalDecomposition().pseudoInverse());
        VectorDimD normal(gbDomain.col(0).cross(gbDomain.col(1)));
        normal.normalize();


        ShiftedDisplacementKernelFT pihatFunction(t,normal);
        auto  lf= LatticeFunction<std::complex<double>,dim-1> (n,basisVectors,pihatFunction);

        // u = f \star b
        for(int i=0; i<dim; ++i) {
            u(i) = bhat[i].dot(lf).real();
            //assert(abs(bhat[i].dot(lf).imag()) < FLT_EPSILON);
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
    DisplacementKernel::DisplacementKernel(const Eigen::Vector<double, Eigen::Dynamic>& _normal):
        Function<DisplacementKernel, double>(),
        normal(_normal)
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

}

#endif