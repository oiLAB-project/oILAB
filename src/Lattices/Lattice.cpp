/* This file is part of gbLAB.
 *
 * gbLAB is distributed without any warranty under the MIT License.
 */


#ifndef gbLAB_Lattice_cpp_
#define gbLAB_Lattice_cpp_

#include <Eigen/Eigenvalues>

#include <LatticeModule.h>
#include <GramMatrix.h>
#include <iomanip>

namespace gbLAB
{


//    /**********************************************************************/
//    template <int dim>
//    typename Lattice<dim>::MatrixDimD Lattice<dim>::getLatticeBasis(const MatrixDimD &A, const MatrixDimD &Q)
//    {
//
//        // Check that Q is orthogonal
//        const MatrixDimD QQT(Q * Q.transpose());
//        const double QQTerror((QQT - MatrixDimD::Identity()).norm());
//        if (QQTerror > 10.0 * DBL_EPSILON * dim * dim)
//        {
//            throw std::runtime_error("The rotation matrix Q is not orthogonal: norm(Q*Q^T-I)="+std::to_string(QQTerror)+"\n");
//        }
//
//        const double Qdet(Q.determinant());
//        if (std::fabs(Q.determinant() - 1.0)>FLT_EPSILON)
//        {
//            throw std::runtime_error("The rotation matrix is not proper: det(Q)="+std::to_string(Qdet)+"\n");
//        }
//        
//        return Q * A;
//    }


    /**********************************************************************/
        /**********************************************************************/
    template <int dim>
    Lattice<dim>::Lattice(const MatrixDimD& A,const MatrixDimD& Fin) :
    /* init */ latticeBasis(Fin*A)
    /* init */,reciprocalBasis(latticeBasis.inverse().transpose())
    /* init */,F(Fin)
    {

    }

    /**********************************************************************/
    template <int dim>
    LatticeDirection<dim> Lattice<dim>::latticeDirection(const VectorDimD &d) const
    {
        const VectorDimD nd(reciprocalBasis.transpose()*d);
        const LatticeVector<dim> temp(LatticeCore<dim>::rationalApproximation(nd),*this);
        const GramMatrix<double,2> G(std::array<VectorDimD,2>{temp.cartesian().normalized(),d.normalized()});
        const double crossNorm(sqrt(G.determinant()));
        if(crossNorm>FLT_EPSILON)
        {
            std::cout<<"input direction="<<d.normalized().transpose()<<std::endl;
            std::cout<<"lattice direction="<<temp.cartesian().normalized().transpose()<<std::endl;
            std::cout<<"cross product norm="<<std::setprecision(15)<<std::scientific<<crossNorm<<std::endl;
            throw std::runtime_error("LATTICE DIRECTION NOT FOUND\n");
        }
        return LatticeDirection<dim>(temp);
    }

    /**********************************************************************/
    template <int dim>
    ReciprocalLatticeDirection<dim> Lattice<dim>::reciprocalLatticeDirection(const VectorDimD &d) const
    {
        const VectorDimD nd(latticeBasis.transpose()*d);
        const ReciprocalLatticeVector<dim> temp(LatticeCore<dim>::rationalApproximation(nd),*this);
        const GramMatrix<double,2> G(std::array<VectorDimD,2>{temp.cartesian().normalized(),d.normalized()});
        const double crossNorm(sqrt(G.determinant()));
        if(crossNorm>FLT_EPSILON)
        {
            std::cout<<"input direction="<<std::setprecision(15)<<std::scientific<<d.normalized().transpose()<<std::endl;
            std::cout<<"reciprocal lattice direction="<<std::setprecision(15)<<std::scientific<<temp.cartesian().normalized().transpose()<<std::endl;
            std::cout<<"cross product norm="<<std::setprecision(15)<<std::scientific<<crossNorm<<std::endl;
            throw std::runtime_error("RECIPROCAL LATTICE DIRECTION NOT FOUND\n");
        }
        return ReciprocalLatticeDirection<dim>(temp);
    }

    /**********************************************************************/
    template <int dim>
    RationalLatticeDirection<dim> Lattice<dim>::rationalLatticeDirection(const VectorDimD &d,
                                                          const typename BestRationalApproximation::LongIntType &maxDen ) const
    {
        const LatticeDirection<dim> ld(latticeDirection(d));
        const BestRationalApproximation bra(d.norm() / ld.cartesian().norm(), maxDen);
        const Rational rat(bra.num, bra.den);
        const RationalLatticeDirection<dim> rld(rat, ld);
        if ((rld.cartesian() - d).squaredNorm() > FLT_EPSILON)
        {
            std::cout << "input vector=" << d.transpose() << std::endl;
            std::cout << "lattice direction=" << ld.cartesian().transpose() << std::endl;
            std::cout << "rational=" << rat << std::endl;
            std::cout << "d.norm()/ld.cartesian().norm()=" << d.norm() / ld.norm() << std::endl;
            throw std::runtime_error("Rational Lattice DirectionType NOT FOUND\n");
        }
        return rld;
    }

    /**********************************************************************/

    template <int dim>
    RationalReciprocalLatticeDirection<dim> Lattice<dim>::rationalReciprocalLatticeDirection(const VectorDimD &d,
                                                                         const typename BestRationalApproximation::LongIntType &maxDen ) const
    {
        const ReciprocalLatticeDirection<dim> rld(reciprocalLatticeDirection(d));
        const BestRationalApproximation bra(d.norm() / rld.cartesian().norm(), maxDen);
        const Rational rat(bra.num, bra.den);
        const RationalReciprocalLatticeDirection<dim> rrld(rat, rld);
        if ((rrld.cartesian() - d).squaredNorm() > FLT_EPSILON)
        {
            std::cout << "input reciprocal vector=" << d.transpose() << std::endl;
            std::cout << "reciprocal lattice direction=" << rld.cartesian().transpose() << std::endl;
            std::cout << "rational=" << rat << std::endl;
            std::cout << "d.norm()/rld.cartesian().norm()=" << d.norm() / rld.norm() << std::endl;
            throw std::runtime_error("Rational Reciprocal Lattice DirectionType NOT FOUND\n");
        }
        return rrld;
    }

    /**********************************************************************/
    template <int dim>
    LatticeVector<dim> Lattice<dim>::latticeVector(const VectorDimD &p) const
    {
        return LatticeVector<dim>(p, *this);
    }

    /**********************************************************************/
    template <int dim>
    ReciprocalLatticeVector<dim> Lattice<dim>::reciprocalLatticeVector(const VectorDimD &p) const
    {
        return ReciprocalLatticeVector<dim>(p, *this);
    }

    /**********************************************************************/
    template<int dim>
    std::vector<LatticeDirection<dim>> Lattice<dim>::planeParallelLatticeBasis(const ReciprocalLatticeDirection<dim>& l,
                                                                               const bool& useRLLL) const
    {
        assert(this == &l.lattice && "Vectors belong to different Lattices.");
        auto outOfPlaneVector = IntegerMath<IntScalarType>::solveBezout(l);
        auto matrix= IntegerMath<IntScalarType>::ccum(outOfPlaneVector);
        std::vector<LatticeDirection<dim>> out;
        int columnIndex= -1;
        for(const auto& column : matrix.colwise())
        {
            columnIndex++;
            if (columnIndex != 0)
                out.push_back(LatticeDirection<dim>(column-column.dot(l)*matrix.col(0),*this));
            else
                out.push_back(LatticeDirection<dim>(column,*this));
        }

        if(!useRLLL) return out;

        Eigen::MatrixXd planeParallelLatticeBasisCartesian(dim,dim-1);
        int index= 0;
        for(auto it=std::next(out.begin()); it!=out.end(); ++it)
        {
            planeParallelLatticeBasisCartesian.col(index)= (*it).cartesian();
            index++;
        }

        planeParallelLatticeBasisCartesian= RLLL(planeParallelLatticeBasisCartesian,0.75).reducedBasis();
        index= 1;
        for(const auto& column : planeParallelLatticeBasisCartesian.colwise())
        {
            out[index]= this->latticeDirection(column);
            index++;
        }
        return out;
    }

    /**********************************************************************/
    template<int dim>
    std::vector<ReciprocalLatticeDirection<dim>> Lattice<dim>::directionOrthogonalReciprocalLatticeBasis(const LatticeDirection<dim>& l,
                                                                                                         const bool& useRLLL) const
    {
        assert(this == &l.lattice && "Vectors belong to different Lattices.");
        auto nonOrthogonalReciprocalVector= IntegerMath<IntScalarType>::solveBezout(l);
        auto matrix= IntegerMath<IntScalarType>::ccum(nonOrthogonalReciprocalVector);
        std::vector<ReciprocalLatticeDirection<dim>> out;
        int columnIndex= -1;
        for(const auto& column : matrix.colwise())
        {
            columnIndex++;
            if (columnIndex != 0) {
                VectorDimI temp= column - column.dot(l) * matrix.col(0);
                out.push_back(ReciprocalLatticeDirection<dim>(ReciprocalLatticeVector<dim>(temp, *this)));
            }
            else {
                VectorDimI temp = column;
                out.push_back(ReciprocalLatticeDirection<dim>(ReciprocalLatticeVector<dim>(temp, *this)));
            }
        }
        if (!useRLLL) return out;

        Eigen::MatrixXd directionOrthogonalReciprocalLatticeBasisCartesian(dim,dim-1);
        int index= 0;
        for(auto it=std::next(out.begin()); it!=out.end(); ++it)
        {
            directionOrthogonalReciprocalLatticeBasisCartesian.col(index)= (*it).cartesian();
            index++;
        }

        directionOrthogonalReciprocalLatticeBasisCartesian= RLLL(directionOrthogonalReciprocalLatticeBasisCartesian,0.75).reducedBasis();
        index= 1;
        for(const auto& column : directionOrthogonalReciprocalLatticeBasisCartesian.colwise())
        {
            out[index]= this->reciprocalLatticeDirection(column);
            index++;
        }
        return out;
    }

    /**********************************************************************/
    template<int dim>
    double Lattice<dim>::interPlanarSpacing(const ReciprocalLatticeDirection<dim>& r) const
    {
        if(&(r.lattice) != this)
            throw(std::runtime_error("The input reciprocal lattice vectors does not belong to the current lattice."));
        return 1.0/r.cartesian().norm();
    }
/*
    template<int dim>
    template<int dm>
    typename std::enable_if<dm==3,std::vector<Lattice<dm>>>::type
    Lattice<dim>::generateCoincidentLattices(const LatticeDirection<dim>& d, const double& maxStrain) const
    {
        std::vector<Lattice<dim>> output;
        return output;
    }
    template<int dim>
    template<int dm>
    typename std::enable_if<dm==2,std::vector<Lattice<dm>>>::type Lattice<dim>::generateCoincidentLattices(const double& maxStrain) const
    {
        std::vector<Lattice<dim>> output;
        return output;
    }
    template<int dim>
    template<int dm>
    typename std::enable_if<(dm!=2) && (dm!=3),std::vector<Lattice<dm>>>::type Lattice<dim>::generateCoincidentLattices(const double& maxStrain) const
    {
        std::vector<Lattice<dim>> output;
        return output;
    }
*/


    template class Lattice<1>;
    template class Lattice<2>;
    template class Lattice<3>;
    template class Lattice<4>;
    template class Lattice<5>;
}
#endif
