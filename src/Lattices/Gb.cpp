//
// Created by Nikhil Chandra Admal on 11/5/22.
//
#ifndef OILAB_GB_CPP
#define OILAB_GB_CPP

#include "Gb.h"

namespace gbLAB
{
    template<int dim>
    Gb<dim>::Gb(const BiCrystal<dim>& bc, const ReciprocalLatticeDirection<dim>& n)
    try :
    /* init */ bc(bc)
    /* init */,nA(&n.lattice == &(bc.A) ? n : bc.getReciprocalLatticeDirectionInA(-1*n.reciprocalLatticeVector()))
    /* init */,nB(&n.lattice == &(bc.B) ? n : bc.getReciprocalLatticeDirectionInB(-1*n.reciprocalLatticeVector()))
    /* init */,basisT(getBasisT(bc,n))
    /* init */,T(Lattice<dim>(bc.dscl.latticeBasis*getBasisT(bc,n).template cast<double>()))
    {
        if (&n.lattice != &(bc.A) && &n.lattice != &(bc.B))
            throw std::runtime_error("The normal does not belong to the reciprocal lattices of A and B  \n");
    }
    catch(std::runtime_error& e)
    {
        std::cout << e.what() << std::endl;
        throw(std::runtime_error("GB construction failed. "));
    }

    template<int dim>
    double Gb<dim>::stepHeightA(const LatticeVector<dim>& d) const
    {
        ReciprocalLatticeDirection<dim> dir= bc.getReciprocalLatticeDirectionInC(nA.reciprocalLatticeVector());
        double cslPlaneSpacing= dir.planeSpacing();
        double step= bc.shiftTensorA(d).cartesian().dot(nA.cartesian().normalized());
        return std::remainder(step,cslPlaneSpacing);
    }

    template<int dim>
    double Gb<dim>::stepHeightB(const LatticeVector<dim>& d) const
    {
        ReciprocalLatticeDirection<dim> dir= bc.getReciprocalLatticeDirectionInC(nA.reciprocalLatticeVector());
        double cslPlaneSpacing= dir.planeSpacing();
        double step= bc.shiftTensorB(d).cartesian().dot(nB.cartesian().normalized());
        return std::remainder(step,cslPlaneSpacing);
    }

    template<int dim> template<int dm>
    typename std::enable_if<dm==2 || dm==3,std::vector<LatticeVector<dim>>>::type
    Gb<dim>::box(std::vector<LatticeVector<dim>>& boxVectors,
                 const double& orthogonality,
                 const int& dsclFactor,
                 std::string filename,
                 bool orient) const
    {
        for (auto iter= std::next(boxVectors.begin()); iter < boxVectors.end(); iter++)
            assert((*iter).dot(bc.getReciprocalLatticeDirectionInC(nA.reciprocalLatticeVector())) == 0 &&
                   "Box vectors not parallel to the grain boundary.");

        auto config= bc.box(boxVectors,orthogonality,dsclFactor);
        std::vector<LatticeVector<dim>> configuration;
        for (const LatticeVector<dim>& latticeVector : config)
        {
            if (&(latticeVector.lattice) == &bc.A && latticeVector.dot(nA)<=0)
                configuration.push_back(latticeVector);
            if (&(latticeVector.lattice) == &bc.B && latticeVector.dot(nB)<=0)
                configuration.push_back(latticeVector);
            if (&(latticeVector.lattice) == &bc.csl || &(latticeVector.lattice) == &bc.dscl)
                configuration.push_back(latticeVector);
        }

        // form the rotation matrix used to orient the system
        MatrixDimD rotation= Eigen::Matrix<double,dim,dim>::Identity();;
        Eigen::Matrix<double,dim,dim-1> orthogonalVectors;
        if (orient) {
            if (dim==3) {
                orthogonalVectors.col(0) = boxVectors[1].cartesian().normalized();
                orthogonalVectors.col(1) = boxVectors[2].cartesian().normalized();
            }
            else if (dim==2)
                orthogonalVectors.col(0)= boxVectors[1].cartesian().normalized();

            rotation=Rotation<dim>(orthogonalVectors);
        }
        assert((rotation*rotation.transpose()).template isApprox(Eigen::Matrix<double,dim,dim>::Identity())
               && "Cannot orient the grain boundary. The GB plane box vectors are not orthogonal.");

        if(!filename.empty()) {
            std::ofstream file;
            file.open(filename);
            if (!file) std::cerr << "Unable to open file";
            file << configuration.size() << std::endl;
            file << "Lattice=\"";

            LatticeVector<dim> origin(-1*boxVectors[0]);
            if (dim == 2) {
                file << (rotation * 2 * boxVectors[0].cartesian()).transpose() << " 0 ";
                file << (rotation * boxVectors[1].cartesian()).transpose() << " 0 ";
                file << " 0 0 1 ";
                file << "\" Properties=atom_types:I:1:pos:R:3:radius:R:1 PBC=\"F T T\" origin=\"";
                file << (rotation * origin.cartesian()).transpose() << " 0.0\"" << std::endl;
                for (const auto &vector: configuration)
                    if (&(vector.lattice) == &bc.A)
                        file << 1 << " " << (rotation * vector.cartesian()).transpose() << " " << 0.0 << "  "
                             << 0.05 << std::endl;
                    else if (&(vector.lattice) == &bc.B)
                        file << 2 << " " << (rotation * vector.cartesian()).transpose() << " " << 0.0 << "  "
                             << 0.05 << std::endl;
                    else if (&(vector.lattice) == &bc.csl)
                        file << 3 << " " << (rotation * vector.cartesian()).transpose() << " " << 0.0 << "  " << 0.2
                             << std::endl;
                    else
                        file << 4 << " " << (rotation * vector.cartesian()).transpose() << " " << 0.0 << "  "
                             << 0.01 << std::endl;
            } else if (dim == 3) {
                file << (rotation * 2 * boxVectors[0].cartesian()).transpose() << " ";
                file << (rotation * boxVectors[1].cartesian()).transpose() << " ";
                file << (rotation * boxVectors[2].cartesian()).transpose() << " ";
                file << "\" Properties=atom_types:I:1:pos:R:3:radius:R:1 PBC=\"F T T\" origin=\"";
                file << (rotation * origin.cartesian()).transpose() << "\"" << std::endl;

                for (const auto &vector: configuration)
                    if (&(vector.lattice) == &bc.A)
                        file << 1 << " " << (rotation * vector.cartesian()).transpose() << "  " << 0.05
                             << std::endl;
                    else if (&(vector.lattice) == &bc.B)
                        file << 2 << " " << (rotation * vector.cartesian()).transpose() << "  " << 0.05
                             << std::endl;
                    else if (&(vector.lattice) == &bc.csl)
                        file << 3 << " " << (rotation * vector.cartesian()).transpose() << "  " << 0.2 << std::endl;
                    else
                        file << 4 << " " << (rotation * vector.cartesian()).transpose() << "  " << 0.01
                             << std::endl;
            }
            file.close();
        }
        return configuration;
    }

    template<int dim> template<int dm>
    typename std::enable_if<dm==2,LatticeVector<dim>>::type
    Gb<dim>::getPeriodVector(const ReciprocalLatticeVector<dim>& axis) const
    {
        LatticeVector<dm> axisAxnA(nA.reciprocalLatticeVector().cross().latticeVector());
        return (LatticeVector<dm>(bc.getLatticeDirectionInC(axisAxnA).latticeVector()));
    }

    template<int dim> template<int dm>
    typename std::enable_if<dm==3,LatticeVector<dim>>::type
    Gb<dim>::getPeriodVector(const ReciprocalLatticeVector<dm>& axis) const
    {
        assert(abs(nA.cartesian().dot(axis.cartesian())) < FLT_EPSILON);
        ReciprocalLatticeDirection<dm> axisA(bc.A.reciprocalLatticeVector(axis.cartesian()));
        LatticeVector<dm> axisAxnA(axisA.reciprocalLatticeVector().cross(nA.reciprocalLatticeVector()).latticeVector());
        return (LatticeVector<dm>(bc.getLatticeDirectionInC(axisAxnA).latticeVector()));
    }


    template<int dim>
    typename Gb<dim>::MatrixDimI Gb<dim>::getBasisT(const BiCrystal<dim> &bc, const ReciprocalLatticeDirection<dim> &n)
    {
        MatrixDimI output;
        MatrixDimI Q = 2*bc.LambdaA-Eigen::Matrix<IntScalarType,dim,dim>::Identity();
        auto nC= bc.getReciprocalLatticeDirectionInC(this->nB.reciprocalLatticeVector());

        MatrixDimI adjM= MatrixDimIExt<IntScalarType,dim>::adjoint(bc.M);
        MatrixDimI adjN= MatrixDimIExt<IntScalarType,dim>::adjoint(bc.N);

        IntScalarType detMN= (bc.M*bc.N).template cast<double>().determinant();
        Eigen::DiagonalMatrix<IntScalarType,dim> nCDiagonal(nC.reciprocalLatticeVector());
        RationalMatrix<dim> temp(adjM*adjN*Q*nCDiagonal, Eigen::Matrix<IntScalarType,dim,dim>::Constant(2*detMN));

        IntScalarType alpha= IntegerMath<IntScalarType>::gcd(temp.integerMatrix.diagonal().eval());
        IntScalarType beta= temp.mu;
        ReciprocalLatticeDirection<dim> m(ReciprocalLatticeVector<dim>(temp.integerMatrix.diagonal().eval(), bc.dscl));

        assert(IntegerMath<IntScalarType>::gcd(alpha,beta) == 1);

        //auto basis= bc.dscl.planeParallelLatticeBasis(m,true);
        auto basis= bc.dscl.planeParallelLatticeBasis(m,false);
        for(int i=0; i<dim; ++i) {
            output.col(i) = basis[i].latticeVector();
            if (i==0) output.col(i)= beta*output.col(i);
        }
        return output;
    }

    template<int dim>
    LatticeVector<dim> Gb<dim>::getLatticeVectorInT(const LatticeVector<dim>& v) const
    {
        MatrixDimI adj= MatrixDimIExt<IntScalarType,dim>::adjoint(basisT);
        IntScalarType det= basisT.template cast<double>().determinant();
        //      basisTA
        // T    -------->   dscl
        if (&(v.lattice) == &(bc.csl)  && IntegerMath<IntScalarType>::gcd(v)%2 == 0)
        {
            VectorDimI integerCoordinates= adj * bc.getLatticeVectorInD(v);
            assert(IntegerMath<IntScalarType>::gcd(integerCoordinates) % det == 0);
            integerCoordinates= integerCoordinates/det;
            return LatticeVector<dim>(integerCoordinates,T);
        }
        else
            throw(std::runtime_error("The input lattice vector should belong to twice the CSL lattice"));

    }

    template<int dim>
    ReciprocalLatticeVector<dim> Gb<dim>::getReciprocalLatticeVectorInT(const ReciprocalLatticeVector<dim>& v) const
    {
        //      basisTA
        // T    -------->   dscl
        if (&(v.lattice) == &(bc.dscl))
            return ReciprocalLatticeVector<dim>((basisT.transpose()*v).eval(),T);
        else
            throw(std::runtime_error("The input reciprocal lattice vector should belong to the reciprocal lattice of the DSCL"));

    }

    template<int dim>
    ReciprocalLatticeDirection<dim> Gb<dim>::getReciprocalLatticeDirectionInT(const ReciprocalLatticeVector<dim>& v) const
    {
        return ReciprocalLatticeDirection<dim>(getReciprocalLatticeVectorInT(v));
    }

    template class Gb<2>;
    template LatticeVector<2> Gb<2>::getPeriodVector<2>(const ReciprocalLatticeVector<2> &axis) const;
    template std::vector<LatticeVector<2>> Gb<2>::box<2>(std::vector<LatticeVector<2>>& boxVectors,
                                                      const double& orthogonality,
                                                      const int& dsclFactor,
                                                      std::string filename,
                                                      bool orient) const;

    template class Gb<3>;
    template LatticeVector<3> Gb<3>::getPeriodVector<3>(const ReciprocalLatticeVector<3> &axis) const;
    template std::vector<LatticeVector<3>> Gb<3>::box<3>(std::vector<LatticeVector<3>>& boxVectors,
                                                         const double& orthogonality,
                                                         const int& dsclFactor,
                                                         std::string filename,
                                                         bool orient) const;

    template class Gb<4>;
    template class Gb<5>;
}


#endif //OILAB_GB_CPP
