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
    LatticeVector<dim> BiCrystal<dim>::getLatticeVectorInA(const LatticeVector<dim> &v) const
    {
        VectorDimI integerCoordinates;

        if(&(v.lattice) == &(this->A))
            return v;
        else if(&(v.lattice) == &(this->csl))
            // U*M*v
            integerCoordinates= this->matrixX() * M * v;
        else
            throw(std::runtime_error("The input lattice vector should belong "
                                     "to lattice A or the CSL"));
        auto temp= LatticeVector<dim>(integerCoordinates,A);
        if (temp.cartesian().dot(v.cartesian()) < 0) temp= -1*temp;
        return temp;
    }

    template<int dim>
    LatticeVector<dim> BiCrystal<dim>::getLatticeVectorInB(const LatticeVector<dim> &v) const
    {
        VectorDimI integerCoordinates;

        if(&(v.lattice) == &(this->B))
            return v;
        else if(&(v.lattice) == &(this->csl))
            // V*N*v
            integerCoordinates= this->matrixV() * N * v;
        else
            throw(std::runtime_error("The input lattice vector should belong "
                                     "to lattice B or the CSL"));
        auto temp= LatticeVector<dim>(integerCoordinates,B);
        if (temp.cartesian().dot(v.cartesian()) < 0) temp= -1*temp;
        return temp;
    }

    template<int dim>
    LatticeVector<dim> BiCrystal<dim>::getLatticeVectorInD(const LatticeVector<dim> &v) const
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
            integerCoordinates = N * M * v;
        else if(&(v.lattice) == &(this->dscl))
            return LatticeVector<dim>(v);
        else
            throw(std::runtime_error("The input lattice vector should belong to one of the four lattices of the bicrystal"));

        auto temp= LatticeVector<dim>(integerCoordinates,dscl);
        if (temp.cartesian().dot(v.cartesian()) < 0) temp= -1*temp;

        //return LatticeVector<dim>(integerCoordinates,dscl);
        return temp;
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

        auto temp= LatticeVector<dim>(integerCoordinates,csl);
        if (temp.cartesian().dot(v.cartesian()) < 0) temp= -1*temp;
        return LatticeDirection<dim>(temp);
    }
    template<int dim>
    LatticeDirection<dim> BiCrystal<dim>::getLatticeDirectionInD(const LatticeVector<dim> &v) const
    {
        return LatticeDirection<dim>(getLatticeVectorInD(v));
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

        auto temp= ReciprocalLatticeVector<dim>(integerCoordinates,A);
        if (temp.cartesian().dot(rv.cartesian()) < 0) temp= -1*temp;
        return ReciprocalLatticeDirection<dim>(temp);
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

        auto temp= ReciprocalLatticeVector<dim>(integerCoordinates,B);
        if (temp.cartesian().dot(rv.cartesian()) < 0) temp= -1*temp;
        return ReciprocalLatticeDirection<dim>(temp);
    }
    template<int dim>
    ReciprocalLatticeDirection<dim> BiCrystal<dim>::getReciprocalLatticeDirectionInC(const ReciprocalLatticeVector<dim>& rv) const
    {
        VectorDimI integerCoordinates;

        if(&(rv.lattice) == &(this->A))
            // M*U^T*rv
            integerCoordinates= M*this->matrixX().transpose()*rv;
        else if(&(rv.lattice) == &(this->B))
            // N*V^T*rv
            integerCoordinates=  N*this->matrixV().transpose()*rv;
        else if(&(rv.lattice) == &(this->csl))
            return ReciprocalLatticeDirection<dim>(rv);
        else if(&(rv.lattice) == &(this->dscl))
            // M*N*rv
            integerCoordinates= M*N*rv;
        else
            throw(std::runtime_error("The input reciprocal lattice vector should belong to one of the four reciprocal lattices of the bicrystal"));

        auto temp= ReciprocalLatticeVector<dim>(integerCoordinates,csl);
        if (temp.cartesian().dot(rv.cartesian()) < 0) temp= -1*temp;
        return ReciprocalLatticeDirection<dim>(temp);
    }
    template<int dim>
    ReciprocalLatticeDirection<dim> BiCrystal<dim>::getReciprocalLatticeDirectionInD(const ReciprocalLatticeVector<dim> &rv) const
    {
        VectorDimI integerCoordinates;

        MatrixDimI adjM= MatrixDimIExt<IntScalarType,dim>::adjoint(M);
        MatrixDimI adjN= MatrixDimIExt<IntScalarType,dim>::adjoint(N);

        if(&(rv.lattice) == &(this->A))
            integerCoordinates= adjN * (this->matrixX().transpose()) * rv;
        else if(&(rv.lattice) == &(this->B))
            integerCoordinates= adjM * (this->matrixV().transpose()) * rv;
        else if(&(rv.lattice) == &(this->csl))
            integerCoordinates= adjN * adjM * rv;
        else
            throw(std::runtime_error("The input reciprocal lattice vector should belong to one of the four reciprocal lattices of the bicrystal"));

        auto temp= ReciprocalLatticeVector<dim>(integerCoordinates,dscl);
        if (temp.cartesian().dot(rv.cartesian()) < 0) temp= -1*temp;
        return ReciprocalLatticeDirection<dim>(temp);
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

    template<int dim> template<int dm>
    typename std::enable_if<dm==2 || dm==3,std::map<typename BiCrystal<dim>::IntScalarType,Gb<dm>>>::type
    BiCrystal<dim>::generateGrainBoundaries(const LatticeDirection<dim>& d, int div) const
    {
        if (&d.lattice != &A && &d.lattice != &B)
            throw std::runtime_error("The tilt axis does not belong to lattices A and B  \n");
        std::vector<Gb<dm>> gbVec;
        double epsilon=1e-8;
        int count= -1;
        IntScalarType keyScale= 1e6;
        auto basis= d.lattice.directionOrthogonalReciprocalLatticeBasis(d,true);
        if (dm==3)
        {
            for (int i = -div; i <= div; ++i)
            {
                for (int j = -div; j <= div; ++j)
                {
                    if (i==0 && j==0) continue;
                    count++;
                    ReciprocalLatticeVector<dm> rv = i * basis[1].reciprocalLatticeVector() + j * basis[2].reciprocalLatticeVector();
                    try
                    {
                        Gb<dm> gb(*this, rv);
                        gbVec.push_back(gb);
                    }
                    catch(std::runtime_error& e)
                    {
                        std::cout << e.what() << std::endl;
                        std::cout << "Unable to form GB with normal = " << rv << std::endl;
                        std::cout << "moving on to next inclination" << std::endl;
                    }
                }
            }
        }
        else if(dm==2)
        {
            auto rv= basis[0].reciprocalLatticeVector();
            gbVec.push_back(Gb<dm>(*this, rv));
        }
        std::map<IntScalarType,Gb<dm>> gbSet;
        for(const Gb<dim>& gb:gbVec)
        {
            double cosAngle;
            cosAngle= gb.nA.cartesian().normalized().dot(gbVec[0].nA.cartesian().normalized());
            if (cosAngle-1>-epsilon) cosAngle= 1.0;
            if (cosAngle+1<epsilon) cosAngle= -1.0;

            double angle= acos(cosAngle);
            IntScalarType key= angle*keyScale;
            gbSet.insert(std::pair<IntScalarType,Gb<dm>>(key,gb));
        }
        return gbSet;
    }

    template<int dim> template<int dm>
    typename std::enable_if<dm==2 || dm==3,std::vector<LatticeVector<dim>>>::type
    BiCrystal<dim>::box(std::vector<LatticeVector<dim>>& boxVectors,
                        const double& orthogonality,
                        const int& dsclFactor,
                        std::string filename,
                        bool orient) const
    {
        assert(orthogonality>=0.0 && orthogonality<=1.0 &&
               "The \"orthogonality\" parameter should be between 0.0 and 1.0");
        assert(dsclFactor>=1 &&
               "The \"dsclFactor\" should be greater than 1.");
        assert(boxVectors.size()==dim);
        for(const auto& boxVector : boxVectors)
        {
            assert(&csl == &boxVector.lattice &&
                   "Box vectors do not belong to the CSL.");
        }

        // Form the box lattice
        MatrixDimD C;
        for (int i=0; i<dim; ++i) {
            C.col(i) = boxVectors[i].cartesian();
        }
        assert(abs(C.determinant()) > FLT_EPSILON && "Box volume is equal to zero.");


        // Adjust boxVector[0] such that it is as orthogonal as possible to boxVector[1]
        auto boxVectorTemp= boxVectors[0];
        ReciprocalLatticeDirection<dim> nC(csl);
        if (dim==2)
            nC= boxVectors[1].cross();
        if (dim==3)
            nC= boxVectors[1].cross(boxVectors[2]);
        auto basis= csl.planeParallelLatticeBasis(nC,true);

        int planesToExplore= nC.stacking();
        MatrixDimI boxLatticeIndices;
        boxLatticeIndices.col(0)= boxVectors[0];
        for (int i=1; i<dim; ++i)
            boxLatticeIndices.col(i)= boxVectors[i]/IntegerMath<long long int>::gcd(boxVectors[i]);
        double minDotProduct= M_PI/2;
        int minStep;

        auto boxVectorUpdated(boxVectors[0]);

        for(int i=0;i<planesToExplore;++i)
        {
            int sign= boxVectors[0].cartesian().dot(basis[0].cartesian())>0? 1 : -1;
            boxVectorTemp= boxVectors[0]+i*sign*basis[0].latticeVector();
            boxLatticeIndices.col(0)= boxVectorTemp;
            Lattice<dim> boxLattice(csl.latticeBasis*boxLatticeIndices.template cast<double>());
            ReciprocalLatticeVector<dim> rC(boxLattice);
            rC(0)=1;
            VectorDimI temp=
                    boxLatticeIndices*boxLattice.planeParallelLatticeBasis(ReciprocalLatticeDirection<dim>(rC),true)[0].latticeVector();
            boxVectorTemp= LatticeVector<dim>(temp,csl);
            double dotProduct= abs(acos(boxVectorTemp.cartesian().normalized().dot(rC.cartesian().normalized())));
            if(dotProduct<minDotProduct) {
                minDotProduct = dotProduct;
                boxVectorUpdated= boxVectorTemp;
                if (dotProduct < (1-orthogonality)*M_PI/2)
                    break;
            }

        }
        boxVectors[0]=boxVectorUpdated;
        C.col(0)= boxVectors[0].cartesian();


        // form the rotation matrix used to orient the system
        MatrixDimD rotation= Eigen::Matrix<double,dim,dim>::Identity();;
        Eigen::Matrix<double,dim,dim-1> orthogonalVectors;
        if (orient) {
            if (dim==3) {
                orthogonalVectors.col(0) = C.col(1).normalized();
                orthogonalVectors.col(1) = C.col(2).normalized();
            }
            else if (dim==2)
                orthogonalVectors.col(0)= C.col(1).normalized();

            rotation=Rotation<dim>(orthogonalVectors);
        }
        assert((rotation*rotation.transpose()).template isApprox(Eigen::Matrix<double,dim,dim>::Identity())
               && "Cannot orient the grain boundary. Box vectors are not orthogonal.");


        std::vector<LatticeVector<dim>> configurationA, configurationB, configurationC, configurationD;
        std::vector<LatticeVector<dim>> configuration;

        std::vector<LatticeVector<dim>> boxVectorsInA, boxVectorsInB, boxVectorsInD;
        // calculate boxVectors in A, B, and D
        for(const auto& boxVector : boxVectors) {
            boxVectorsInA.push_back(getLatticeVectorInA(boxVector));
            boxVectorsInB.push_back(getLatticeVectorInB(boxVector));
            boxVectorsInD.push_back(getLatticeVectorInD(boxVector));
        }

        // prepare boxVectors for D
        auto dsclVector=getLatticeDirectionInD(boxVectors[0]).latticeVector();
        auto nD= getReciprocalLatticeDirectionInD(nC.reciprocalLatticeVector());
        if(abs((dsclFactor*dsclVector).dot(nD)) < abs(boxVectorsInD[0].dot(nD)))
            boxVectorsInD[0]= dsclFactor*dsclVector;

        std::vector<LatticeVector<dim>> boxVectorsForA(boxVectorsInA),
                boxVectorsForB(boxVectorsInB),
                boxVectorsForC(boxVectors),
                boxVectorsForD(boxVectorsInD);
        boxVectorsForA[0]=2*boxVectorsInA[0];
        boxVectorsForB[0]=2*boxVectorsInB[0];
        boxVectorsForC[0]=2*boxVectors[0];
        boxVectorsForD[0]=2*boxVectorsInD[0];

        configurationA= A.box(boxVectorsForA);
        configurationB= B.box(boxVectorsForB);
        configurationC= csl.box(boxVectorsForC);
        configurationD= dscl.box(boxVectorsForD);

        LatticeVector<dim> origin(-1*boxVectors[0]);
        for(auto& vector : configurationA)
            vector= vector + LatticeVector<dim>(-1*boxVectorsInA[0]);
        for(auto& vector : configurationB)
            vector= vector + LatticeVector<dim>(-1*boxVectorsInB[0]);
        for(auto& vector : configurationC)
            vector= vector + origin;
        for(auto& vector : configurationD)
            vector= vector + LatticeVector<dim>(-1*boxVectorsInD[0]);

        configuration= configurationA;
        configuration.insert(configuration.end(),configurationB.begin(),configurationB.end());
        configuration.insert(configuration.end(),configurationC.begin(),configurationC.end());
        configuration.insert(configuration.end(),configurationD.begin(),configurationD.end());

        if(!filename.empty()) {
            std::ofstream file;
            file.open(filename);
            if (!file) std::cerr << "Unable to open file";
            file << configuration.size() << std::endl;
            file << "Lattice=\"";

            if (dim==2) {
                file << (rotation*boxVectorsForC[0].cartesian()).transpose() << " 0 ";
                file << (rotation*boxVectorsForC[1].cartesian()).transpose() << " 0 ";
                file << " 0 0 1 ";
                file << "\" Properties=atom_types:I:1:pos:R:3:radius:R:1 PBC=\"F T T\" origin=\"";
                file << (rotation*origin.cartesian()).transpose() << " 0.0\"" << std::endl;
                for (const auto &vector: configurationA)
                    file << 1 << " " << (rotation*vector.cartesian()).transpose() << " " << 0.0 << "  " << 0.05 << std::endl;
                for (const auto &vector: configurationB)
                    file << 2 << " " << (rotation*vector.cartesian()).transpose() << " " << 0.0 << "  " << 0.05 << std::endl;
                for (const auto &vector: configurationC)
                    file << 3 << " " << (rotation*vector.cartesian()).transpose() << " " << 0.0 << "  " << 0.2 << std::endl;
                for (const auto &vector: configurationD)
                    file << 4 << " " << (rotation*vector.cartesian()).transpose() << " " << 0.0 << "  " << 0.01 << std::endl;
            }
            else if (dim==3){
                file << (rotation*boxVectorsForC[0].cartesian()).transpose()  << " ";
                file << (rotation*boxVectorsForC[1].cartesian()).transpose() << " ";
                file << (rotation*boxVectorsForC[2].cartesian()).transpose() << " ";
                file << "\" Properties=atom_types:I:1:pos:R:3:radius:R:1 PBC=\"F T T\" origin=\"";
                file << (rotation*origin.cartesian()).transpose() << "\"" << std::endl;

                for (const auto &vector: configurationA)
                    file << 1 << " " << (rotation*vector.cartesian()).transpose() << "  " << 0.05 << std::endl;
                for (const auto &vector: configurationB)
                    file << 2 << " " << (rotation*vector.cartesian()).transpose() << "  " << 0.05 << std::endl;
                for (const auto &vector: configurationC)
                    file << 3 << " " << (rotation*vector.cartesian()).transpose() << "  " << 0.2 << std::endl;
                for (const auto &vector: configurationD)
                    file << 4 << " " << (rotation*vector.cartesian()).transpose() << "  " << 0.01 << std::endl;
            }


            file.close();
        }
        return configuration;
    }

//    template class BiCrystal<1>;
    template class BiCrystal<2>;
    template std::map<BiCrystal<2>::IntScalarType, Gb<2>>
        BiCrystal<2>::generateGrainBoundaries<2>(const LatticeDirection<2> &d, int div) const;
    template std::vector<LatticeVector<2>>
            BiCrystal<2>::box<2>(std::vector<LatticeVector<2>> &boxVectors,
                                 const double &orthogonality, const int &dsclFactor,
                                 std::string filename, bool orient) const;

    template class BiCrystal<3>;
    template std::map<BiCrystal<3>::IntScalarType, Gb<3>>
        BiCrystal<3>::generateGrainBoundaries<3>(const LatticeDirection<3> &d, int div) const;
    template std::vector<LatticeVector<3>>
    BiCrystal<3>::box<3>(std::vector<LatticeVector<3>> &boxVectors,
                         const double &orthogonality, const int &dsclFactor,
                         std::string filename, bool orient) const;

    template class BiCrystal<4>;
    template class BiCrystal<5>;

} // end namespace
#endif

