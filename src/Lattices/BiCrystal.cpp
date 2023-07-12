/* This file is part of gbLAB.
 *
 * gbLAB is distributed without any warranty under the MIT License.
 */


#ifndef gbLAB_BiCrystal_cpp_
#define gbLAB_BiCrystal_cpp_

#include <LatticeModule.h>

namespace gbLAB
{
        
    template <int dim>
    typename BiCrystal<dim>::MatrixDimI BiCrystal<dim>::getM(const RationalMatrix<dim>& rm,
                        const SmithDecomposition<dim>& sd)
    {
        typename BiCrystal<dim>::MatrixDimI M=BiCrystal<dim>::MatrixDimI::Identity();
        for(int i=0;i<dim;++i)
        {
            const auto& dii(sd.matrixD()(i,i));
            M(i,i)=dii/IntegerMath<IntScalarType>::gcd(rm.mu,dii);
        }
        return M;
    }

    template <int dim>
    typename BiCrystal<dim>::MatrixDimI BiCrystal<dim>::getN(const RationalMatrix<dim>& rm,
                        const SmithDecomposition<dim>& sd)
    {
        typename BiCrystal<dim>::MatrixDimI N=BiCrystal<dim>::MatrixDimI::Identity();
        for(int i=0;i<dim;++i)
        {
            const auto& dii=sd.matrixD()(i,i);
            N(i,i)=rm.mu/IntegerMath<IntScalarType>::gcd(rm.mu,dii);
        }
        return N;
    }

    template <int dim>
    typename BiCrystal<dim>::MatrixDimI BiCrystal<dim>::getLambdaA(const typename BiCrystal<dim>::MatrixDimI& M,
                                                                   const typename BiCrystal<dim>::MatrixDimI& N)
    {
        BiCrystal<dim>::MatrixDimI LambdaA;
        for(int col=0; col<dim; ++col)
        {
            BiCrystal<dim>::VectorDimI x,y;
            for (int i = 0; i < dim; ++i)
            {
                IntScalarType a = N(i, i);
                IntScalarType b = -M(i, i);
                IntScalarType c = -(i==col);
                IntegerMath<IntScalarType>::solveDiophantine2vars(a, b, c, x(i), y(i));
            }
            LambdaA.col(col)= M*y;
        }
        return LambdaA;
    }

    template <int dim>
    typename BiCrystal<dim>::MatrixDimI BiCrystal<dim>::getLambdaB(const typename BiCrystal<dim>::MatrixDimI& M,
                                                                   const typename BiCrystal<dim>::MatrixDimI& N)
    {
        typename BiCrystal<dim>::MatrixDimI LambdaB;
        for(int col=0; col<dim; ++col)
        {
            BiCrystal<dim>::VectorDimI x,y;
            for (int i = 0; i < dim; ++i)
            {
                IntScalarType a = N(i, i);
                IntScalarType b = -M(i, i);
                IntScalarType c = (i==col);
                IntegerMath<IntScalarType>::solveDiophantine2vars(a, b, c, x(i), y(i));
            }
            LambdaB.col(col)= N*x;
        }
        return LambdaB;
    }

    template <int dim>
    typename BiCrystal<dim>::MatrixDimD BiCrystal<dim>::getCSLBasis(const Lattice<dim>& A,
                                                                             const Lattice<dim>& B,
                                                                             const SmithDecomposition<dim>& sd,
                                                                             const typename BiCrystal<dim>::MatrixDimI& M,
                                                                             const typename BiCrystal<dim>::MatrixDimI& N,
                                                                             const bool& useRLLL)
    {
        // The transition matrix is T=P/sigma, where P=rm.integerMatrix is
        // an integer matrix and sigma=rm.sigma is an integer
        // The integer matrix P can be decomposed as P=X*D*Y using the Smith decomposition.
        // X and Y are unimodular, and D is diagonal with D(k,k) dividing D(k+1,k+1)
        // The decomposition also computes the matices U and V such that D=U*P*V
        // From T=inv(A)*B=P/sigma=X*D*Y/sigma=X*D*inv(V)/sigma, we have
        // B1*(sigma*I)=A1*D
        // where
        // B1=B*V
        // A1=A*X
        // Since V and X are unimodular matrices, B1 and A1 are new bases
        // of the lattices B and A, respectively. Moreover, since
        // (sigma*I) and D are diagonal, the columns of B1 and A1 are
        // proportional, with rational proportionality factors different for each column.
        // For example, columns "i" read
        // b1_i*sigma=a1_i*D(i,i)
        // Therefore, the i-th primitive vectors of the CSL is
        // c_i=b1_i*sigma/gcd(sigma,D(i,i))=a1_i*D(i,i)/gcd(sigma,D(i,i))
        // or, in matrix form
        // C=B1*N=A1*M, that is
        // C=B*V*N=A*X*M
        // where M=diag(D(i,i)/gcd(sigma,D(i,i))) and
        //       N=diag(sigma/gcd(sigma,D(i,i))) and

        const auto C1(A.latticeBasis*(sd.matrixX()*M).template cast<double>());
        const auto C2(B.latticeBasis*(sd.matrixV()*N).template cast<double>());
        if ((C1-C2).norm()/C1.norm()>FLT_EPSILON || (C1-C2).norm()/C2.norm()>FLT_EPSILON)
        {
            throw std::runtime_error("CSL calculation failed.\n");
        }

        if(useRLLL)
        {
            return RLLL(0.5*(C1+C2),0.75).reducedBasis();
        }
        else
        {
            return 0.5*(C1+C2);
        }
    }

    template <int dim>
    typename BiCrystal<dim>::MatrixDimD BiCrystal<dim>::getDSCLBasis(const Lattice<dim>& A,
                                   const Lattice<dim>& B,
                                   const SmithDecomposition<dim>& sd,
                                   const typename BiCrystal<dim>::MatrixDimI& M,
                                   const typename BiCrystal<dim>::MatrixDimI& N,
                                   const bool& useRLLL)
    {

        const auto D1(A.latticeBasis*sd.matrixX().template cast<double>()*N.template cast<double>().inverse());
        const auto D2(B.latticeBasis*sd.matrixV().template cast<double>()*M.template cast<double>().inverse());
        if ((D1-D2).norm()/D1.norm()>FLT_EPSILON || (D1-D2).norm()/D2.norm()>FLT_EPSILON)
        {
            throw std::runtime_error("DSCL calculation failed.\n");
        }
        if(useRLLL)
        {
            return RLLL(0.5*(D1+D2),0.75).reducedBasis();
        }
        else
        {
            return 0.5*(D1+D2);
        }
    }


    template <int dim>
    BiCrystal<dim>::BiCrystal(const Lattice<dim>& A_in,
              const Lattice<dim>& B_in,
              const bool& useRLLL) try :
    /* init */ RationalMatrix<dim>(A_in.reciprocalBasis.transpose()*B_in.latticeBasis)
    /* init */,SmithDecomposition<dim>(this->integerMatrix)
    /* init */,A(A_in)
    /* init */,B(B_in)
    /* init */,M(getM(*this,*this))
    /* init */,N(getN(*this,*this))
    /* init */,sigmaA(round(M.template cast<double>().determinant()))
    /* init */,sigmaB(round(N.template cast<double>().determinant()))
    /* init */,sigma(std::abs(sigmaA)==std::abs(sigmaB)? std::abs(sigmaA) : 0)
    /* init */, csl(getCSLBasis (A,B,*this,M,N,useRLLL),MatrixDimD::Identity())
    /* init */,dscl(getDSCLBasis(A,B,*this,M,N,useRLLL),MatrixDimD::Identity())
    /* init */,Ap(A.latticeBasis*this->matrixX().template cast<double>())
    /* init */,Bp(B.latticeBasis*this->matrixV().template cast<double>())
    /* init */,LambdaA(getLambdaA(M,N))
    /* init */,LambdaB(getLambdaB(M,N))
    {

        if(true)
        {//verify that CSL can be obtained as multiple of A and B

            const MatrixDimD tempA(A.reciprocalBasis.transpose()*csl.latticeBasis);
            if ((tempA-tempA.array().round().matrix()).norm()/tempA.norm()>FLT_EPSILON)
            {
                //std::cout << (tempA-tempA.array().round().matrix()).norm()/tempA.norm() << std::endl;
                throw std::runtime_error("CSL is not a multiple of lattice A.\n");
            }

            const MatrixDimD tempB(B.reciprocalBasis.transpose()*csl.latticeBasis);
            if ((tempB-tempB.array().round().matrix()).norm()/tempB.norm()>FLT_EPSILON)
            {
                //std::cout << (tempB-tempB.array().round().matrix()).norm()/tempB.norm() << std::endl;
                throw std::runtime_error("CSL is not a multiple of lattice B\n");
            }
        }

        if(true)
        {//verify that A and B are multiples of DSCL

            const MatrixDimD tempA(dscl.reciprocalBasis.transpose()*A.latticeBasis);
            if ((tempA-tempA.array().round().matrix()).norm()/tempA.norm()>FLT_EPSILON)
            {
                //std::cout << (tempA-tempA.array().round().matrix()).norm()/tempA.norm() << std::endl;
                throw std::runtime_error("Lattice A is not a multiple of the DSCL\n");
            }

            const MatrixDimD tempB(dscl.reciprocalBasis.transpose()*B.latticeBasis);
            if ((tempB-tempB.array().round().matrix()).norm()/tempB.norm()>FLT_EPSILON)
            {
                //std::cout << (tempB-tempB.array().round().matrix()).norm()/tempB.norm() << std::endl;
                throw std::runtime_error("Lattice B is not a multiple of the DSCL\n");
            }
        }

        if(true)
        {//verify LambdaA + LambdaB = I
            if (!(LambdaA+LambdaB).isIdentity())
            {
                throw std::runtime_error("LambdaA + LambdaB != I\n");
            }
        }

    }

    catch(std::runtime_error& e)
    {
        std::cout << e.what() << std::endl;
        throw(std::runtime_error("Bicrystal construction failed. "));
    }

    template<int dim>
    LatticeDirection<dim> BiCrystal<dim>::getLatticeDirectionInC(const LatticeVector<dim> &v) const
    {
        VectorDimI integerCoordinates;

        MatrixDimI adjX= MatrixDimIExt<IntScalarType,dim>::adjoint(this->matrixX());
        MatrixDimI adjV= MatrixDimIExt<IntScalarType,dim>::adjoint(this->matrixV());
        MatrixDimI adjM= MatrixDimIExt<IntScalarType,dim>::adjoint(M);
        MatrixDimI adjN= MatrixDimIExt<IntScalarType,dim>::adjoint(N);
        MatrixDimI adjMN= MatrixDimIExt<IntScalarType,dim>::adjoint(M*N);

        if(&(v.lattice) == &(this->A))
            // inv(M)*inv(U)*v
            integerCoordinates = adjM * adjX * v;
        else if(&(v.lattice) == &(this->B))
            // inv(N)*inv(V)*v
            integerCoordinates = adjN * adjV * v;
        else if(&(v.lattice) == &(this->csl))
            return LatticeDirection<dim>(v);
        else if(&(v.lattice) == &(this->dscl))
            integerCoordinates = adjMN* v;
        else
            throw(std::runtime_error("The input reciprocal lattice vector should belong to one of the four reciprocal lattices of the bicrystal"));

        return LatticeDirection<dim>(LatticeVector<dim>(integerCoordinates,csl));
    }
    template<int dim>
    LatticeDirection<dim> BiCrystal<dim>::getLatticeDirectionInD(const LatticeVector<dim> &v) const
    {
        VectorDimI integerCoordinates;

        MatrixDimI adjX= MatrixDimIExt<IntScalarType,dim>::adjoint(this->matrixX());
        MatrixDimI adjV= MatrixDimIExt<IntScalarType,dim>::adjoint(this->matrixV());

        if(&(v.lattice) == &(this->A))
            // N*inv(U)*v
            integerCoordinates = N * adjX * v;
        else if(&(v.lattice) == &(this->B))
            // M*inv(V)*v
            integerCoordinates = M * adjV * v;
        else if(&(v.lattice) == &(this->csl))
            // N*M*v
            return LatticeDirection<dim>(LatticeVector<dim>((VectorDimI) (N*M*v),dscl));
        else if(&(v.lattice) == &(this->dscl))
            return LatticeDirection<dim>(v);
        else
            throw(std::runtime_error("The input reciprocal lattice vector should belong to one of the four reciprocal lattices of the bicrystal"));

        return LatticeDirection<dim>(LatticeVector<dim>(integerCoordinates,dscl));
    }

    template<int dim>
    ReciprocalLatticeDirection<dim> BiCrystal<dim>::getReciprocalLatticeDirectionInA(const ReciprocalLatticeVector<dim>& rv) const
    {
        VectorDimI integerCoordinates;

        MatrixDimI adjX = MatrixDimIExt<IntScalarType, dim>::adjoint(this->matrixX());
        MatrixDimI adjM= MatrixDimIExt<IntScalarType,dim>::adjoint(M);

        if(&(rv.lattice) == &(this->A))
            return ReciprocalLatticeDirection<dim>(rv);
        else if(&(rv.lattice) == &(this->B))
            // U^-T*inverse(M)*N*V^T
            integerCoordinates= adjX.transpose() * adjM * N * (this->matrixV()).transpose() * rv;
        else if(&(rv.lattice) == &(this->csl))
            // U^-T*inverse(M)
            integerCoordinates= adjX.transpose() * adjM * rv;
        else if(&(rv.lattice) == &(this->dscl))
            // U^-T*N*rv
            integerCoordinates = adjX.transpose() * N * rv;
        else
            throw(std::runtime_error("The input reciprocal lattice vector should belong to one of the four reciprocal lattices of the bicrystal"));
        return ReciprocalLatticeDirection<dim>(ReciprocalLatticeVector<dim>(integerCoordinates,A));
    }
    template<int dim>
    ReciprocalLatticeDirection<dim> BiCrystal<dim>::getReciprocalLatticeDirectionInB(const ReciprocalLatticeVector<dim>& rv) const
    {
        VectorDimI integerCoordinates;

        MatrixDimI adjX= MatrixDimIExt<IntScalarType,dim>::adjoint(this->matrixX());
        MatrixDimI adjV= MatrixDimIExt<IntScalarType,dim>::adjoint(this->matrixV());
        MatrixDimI adjM= MatrixDimIExt<IntScalarType,dim>::adjoint(M);
        MatrixDimI adjN= MatrixDimIExt<IntScalarType,dim>::adjoint(N);

        if(&(rv.lattice) == &(this->A))
            // V^-T*inverse(N)*M*U^T
            integerCoordinates= adjV.transpose() * adjN * M * (this->matrixX()).transpose() * rv;
        else if(&(rv.lattice) == &(this->B))
            return ReciprocalLatticeDirection<dim>(rv);
        else if(&(rv.lattice) == &(this->csl))
            // V^-T*inverse(N)*rv
            integerCoordinates= adjV.transpose() * adjN * rv;
        else if(&(rv.lattice) == &(this->dscl))
            // V^-T*M*rv
            integerCoordinates= adjV.transpose() * M * rv;
        else
            throw(std::runtime_error("The input reciprocal lattice vector should belong to one of the four reciprocal lattices of the bicrystal"));
        return ReciprocalLatticeDirection<dim>(ReciprocalLatticeVector<dim>(integerCoordinates,B));
    }
    template<int dim>
    ReciprocalLatticeDirection<dim> BiCrystal<dim>::getReciprocalLatticeDirectionInC(const ReciprocalLatticeVector<dim>& rv) const
    {
        if(&(rv.lattice) == &(this->A))
            // M*U^T*rv
            return ReciprocalLatticeDirection<dim>(ReciprocalLatticeVector<dim>((VectorDimI) (M*this->matrixX().transpose()*rv),csl));
        else if(&(rv.lattice) == &(this->B))
            // N*V^T*rv
            return ReciprocalLatticeDirection<dim>(ReciprocalLatticeVector<dim>((VectorDimI) (N*this->matrixV().transpose()*rv),csl));
        else if(&(rv.lattice) == &(this->csl))
            return ReciprocalLatticeDirection<dim>(rv);
        else if(&(rv.lattice) == &(this->dscl))
            // M*N*rv
            return ReciprocalLatticeDirection<dim>(ReciprocalLatticeVector<dim>((VectorDimI) (M*N*rv),csl));
        else
            throw(std::runtime_error("The input reciprocal lattice vector should belong to one of the four reciprocal lattices of the bicrystal"));
    }
    template<int dim>
    ReciprocalLatticeDirection<dim> BiCrystal<dim>::getReciprocalLatticeDirectionInD(const ReciprocalLatticeVector<dim> &rv) const
    {
        VectorDimI integerCoordinates;

        MatrixDimI adjM= MatrixDimIExt<IntScalarType,dim>::adjoint(M);
        MatrixDimI adjN= MatrixDimIExt<IntScalarType,dim>::adjoint(N);

        if(&(rv.lattice) == &(this->A))
            integerCoordinates= adjN * rv;
        else if(&(rv.lattice) == &(this->B))
            integerCoordinates= adjM * rv;
        else
            throw(std::runtime_error("The input reciprocal lattice vector should belong to one of the four reciprocal lattices of the bicrystal"));

        return ReciprocalLatticeDirection<dim>(ReciprocalLatticeVector<dim>(integerCoordinates,dscl));
    }

    template<int dim>
    LatticeVector<dim> BiCrystal<dim>::shiftTensorA(const LatticeVector<dim>& d) const
    {
        if(&d.lattice != &this->dscl)
            throw(std::runtime_error("Input vector is not a DSCL vectors"));
        return LatticeVector<dim>((LambdaA*d).eval(),d.lattice);
    }

    template<int dim>
    LatticeVector<dim> BiCrystal<dim>::shiftTensorB(const LatticeVector<dim>& d) const
    {
        if(&d.lattice != &this->dscl)
            throw(std::runtime_error("Input vector is not a DSCL vectors"));
        return LatticeVector<dim>((LambdaB*d).eval(),d.lattice);
    }


//    template class BiCrystal<1>;
    template class BiCrystal<2>;
    template class BiCrystal<3>;
    template class BiCrystal<4>;
    template class BiCrystal<5>;

} // end namespace
#endif

