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
        RLLL rlll(latticeBasis,0.75);
        auto structureMatrix= rlll.reducedBasis();
        auto U= rlll.unimodularMatrix();
        Lattice<dim> reducedLattice(structureMatrix);
        VectorDimD nd;
        nd= reducedLattice.reciprocalBasis.transpose()*d;
        LatticeVector<dim> tempReduced(LatticeCore<dim>::rationalApproximation(nd),reducedLattice);
        VectorDimI tempRecovered= U*tempReduced;
        LatticeVector<dim> temp(tempRecovered, *this);

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
        RLLL rlll(latticeBasis,0.75);
        auto structureMatrix= rlll.reducedBasis();
        MatrixDimI U= rlll.unimodularMatrix();
        Lattice<dim> reducedLattice(structureMatrix);
        VectorDimD nd;
        nd= reducedLattice.latticeBasis.transpose()*d;
        ReciprocalLatticeVector<dim> tempReduced(LatticeCore<dim>::rationalApproximation(nd),reducedLattice);
        VectorDimI tempRecovered(MatrixDimIExt<IntScalarType,dim>::adjoint(U.matrix()).transpose()* tempReduced);
        ReciprocalLatticeVector<dim> temp(tempRecovered, *this);

        const GramMatrix<double,2> G(std::array<VectorDimD,2>{temp.cartesian().normalized(),d.normalized()});
        const double crossNorm(sqrt(abs(G.determinant())));
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
            std::cout << "d.norm()/ld.cartesian().norm()=" << d.norm() / ld.latticeVector().norm() << std::endl;
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
            std::cout << "d.norm()/rld.cartesian().norm()=" << d.norm() / rld.reciprocalLatticeVector().norm() << std::endl;
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
        auto outOfPlaneVector = IntegerMath<IntScalarType>::solveBezout(l.reciprocalLatticeVector());
        auto matrix= IntegerMath<IntScalarType>::ccum(outOfPlaneVector);
        std::vector<LatticeDirection<dim>> out;
        int columnIndex= -1;
        for(const auto& column : matrix.colwise())
        {
            columnIndex++;
            if (columnIndex != 0)
                out.push_back(LatticeDirection<dim>(column-column.dot(l.reciprocalLatticeVector())*matrix.col(0),*this));
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

        if(dim>2) {
            planeParallelLatticeBasisCartesian = RLLL(planeParallelLatticeBasisCartesian, 0.75).reducedBasis();
            index = 1;
            for (const auto &column: planeParallelLatticeBasisCartesian.colwise()) {
                out[index] = this->latticeDirection(column);
                index++;
            }
        }

        // (A^T A)^{-1} A^T
        Eigen::MatrixXd pseudoInverse(dim-1,dim);
        pseudoInverse= (planeParallelLatticeBasisCartesian.transpose() * planeParallelLatticeBasisCartesian).inverse() *
                       (planeParallelLatticeBasisCartesian.transpose());

        LatticeVector<dim> temp(*this);
        for(int i=0;i<dim-1;i++)
        {
            temp= temp + round(out[0].cartesian().dot(pseudoInverse.row(i)))*out[i+1].latticeVector();
        }
        out[0]= LatticeDirection<dim>(out[0].latticeVector()-temp);


        return out;
    }

    /**********************************************************************/
    template<int dim>
    std::vector<ReciprocalLatticeDirection<dim>> Lattice<dim>::directionOrthogonalReciprocalLatticeBasis(const LatticeDirection<dim>& l,
                                                                                                         const bool& useRLLL) const
    {
        assert(this == &l.lattice && "Vectors belong to different Lattices.");
        auto nonOrthogonalReciprocalVector= IntegerMath<IntScalarType>::solveBezout(l.latticeVector());
        auto matrix= IntegerMath<IntScalarType>::ccum(nonOrthogonalReciprocalVector);
        std::vector<ReciprocalLatticeDirection<dim>> out;
        int columnIndex= -1;
        for(const auto& column : matrix.colwise())
        {
            columnIndex++;
            if (columnIndex != 0) {
                VectorDimI temp= column - column.dot(l.latticeVector()) * matrix.col(0);
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

    template<int dim> template<int dm>
    typename std::enable_if<dm==3,std::vector<typename Lattice<dim>::MatrixDimD>>::type
    Lattice<dim>::generateCoincidentLattices(const ReciprocalLatticeDirection<dim>& rd, const double& maxDen, const int& N) const
    {
        std::vector<MatrixDimD> output;
        std::map<IntScalarType,MatrixDimD> temp;
        auto basis= planeParallelLatticeBasis(rd);
        double epsilon=1e-8;
        IntScalarType keyScale= 1e6;

        auto b1= basis[1].cartesian();
        auto b2= basis[2].cartesian();

        // Following the notation in Algorithm 3 of
        // "Interface dislocations and grain boundary disconnections using Smith normal bicrystallography"
        auto q1= b1;
        std::vector<std::pair<int,int>> coPrimePairs= farey(N,false);
        for(const auto& pair : coPrimePairs)
        {
            VectorDimI pIndices;
            auto q2= pair.first*b1+pair.second*b2;

            if (abs(q2.norm())<epsilon ) continue;
            double ratio= q2.norm()/q1.norm();
            BestRationalApproximation bra(ratio,maxDen);
            double error= ratio-static_cast<double>(bra.num)/bra.den;
            if (abs(error) > epsilon)
                continue;
            else
            {
                double cosTheta= b1.dot(q2)/(q1.norm()*q2.norm());
                if (cosTheta-1>-epsilon) cosTheta= 1.0;
                if (cosTheta+1<epsilon) cosTheta= -1.0;
                double theta= acos(cosTheta);
                MatrixDimD rotation;
                rotation= Eigen::AngleAxis<double>(theta,rd.cartesian().normalized());
                IntScalarType key= theta*keyScale;
                temp.insert(std::pair<IntScalarType,MatrixDimD>(key,rotation));
            }
        }
        std::transform(temp.begin(), temp.end(),
                       std::back_inserter(output),
                       [](const std::pair<IntScalarType,MatrixDimD>& p) {
                           return p.second;
                       });
        return output;
    }


    template<int dim> template<int dm>
    typename std::enable_if<dm==2,std::vector<typename Lattice<dim>::MatrixDimD>>::type
    Lattice<dim>::generateCoincidentLattices(const double& maxStrain, const int& maxDen, const int& N) const
    {
        std::vector<MatrixDimD> output(generateCoincidentLattices(*this,maxStrain,maxDen,N));
        return output;
    }


    template<int dim> template<int dm>
    typename std::enable_if<dm==2,std::vector<typename Lattice<dim>::MatrixDimD>>::type
    Lattice<dim>::generateCoincidentLattices(const Lattice<dim>& undeformedLattice,
                                             const double& maxStrain,
                                             const int& maxDen,
                                             const int& N) const
    {
        int numberOfConfigurations= 0;
        const int maxConfigurations= 80;
        std::vector<MatrixDimD> output;
        std::map<IntScalarType,MatrixDimD> temp;
        MatrixDimI mn(MatrixDimI::Zero());
        MatrixDimI md(MatrixDimI::Ones());
        std::vector<std::pair<int,int>> coPrimePairs= farey(N,false);
        double epsilon=1e-8;
        IntScalarType keyScale= 1e6;

        for(const auto& pair1 : coPrimePairs)
        {
            VectorDimI pIndices;
            pIndices << pair1.first, pair1.second;
            LatticeVector<dm> q2(pIndices,undeformedLattice);
            double ratio= q2.cartesian().norm() / latticeBasis.col(0).norm();


            if (abs(maxStrain) < epsilon)
            {
                BestRationalApproximation alpha(ratio, maxDen);
                double error = ratio - static_cast<double>(alpha.num) / alpha.den;
                if (abs(error) > epsilon)
                    continue;
                else
                {
                    double cosTheta = latticeBasis.col(0).dot(q2.cartesian()) / (latticeBasis.col(0).norm() * q2.cartesian().norm());
                    if (cosTheta - 1 > -epsilon) cosTheta = 1.0;
                    if (cosTheta + 1 < epsilon) cosTheta = -1.0;
                    double theta = acos(cosTheta);
                    Eigen::Rotation2D<double> rotation(theta);
                    IntScalarType key= theta*keyScale;
                    temp.insert(std::pair<IntScalarType,MatrixDimD>(key,rotation.toRotationMatrix()));
                }
            }
            else
            {
                RationalApproximations <IntScalarType> alphaSequence(ratio, maxDen, ratio*maxStrain);
                for (const auto& alpha: alphaSequence.approximations)
                {
                    RationalLatticeDirection<dm> q2ByAlpha(Rational<IntScalarType>(alpha.d, alpha.n), q2);

                    mn.col(0) = q2 * alpha.d;
                    md.col(0).setConstant(alpha.n);

                    double s1 =
                            (q2ByAlpha.cartesian().squaredNorm() - latticeBasis.col(0).squaredNorm()) /
                            latticeBasis.col(0).squaredNorm();
                    if (abs(s1) > maxStrain) continue;

                    for (const auto &pair2: coPrimePairs) {
                        VectorDimI qIndices;
                        qIndices << pair2.first, pair2.second;
                        LatticeVector<dm> r2(qIndices, undeformedLattice);
                        double ratio2= r2.cartesian().norm() / latticeBasis.col(1).norm();
                        RationalApproximations<IntScalarType> betaSequence(ratio2, maxDen,maxStrain*ratio2);
                        for(const auto& beta : betaSequence.approximations) {
                            RationalLatticeDirection<dm> r2ByBeta(Rational<IntScalarType>(beta.d, beta.n), r2);
                            double s2 = (r2ByBeta.cartesian().squaredNorm() - latticeBasis.col(1).squaredNorm()) /
                                        latticeBasis.col(1).squaredNorm();
                            double s3 = (q2ByAlpha.cartesian().dot(r2ByBeta.cartesian()) -
                                         latticeBasis.col(0).dot(latticeBasis.col(1))) /
                                        (latticeBasis.col(0).norm() * latticeBasis.col(1).norm());
                            if (abs(s2) > maxStrain || abs(s3) > maxStrain) continue;

                            // calculate the deformation gradient
                            MatrixDimD F = q2ByAlpha.cartesian() * reciprocalBasis.col(0).transpose() +
                                           r2ByBeta.cartesian() * reciprocalBasis.col(1).transpose();
                            if (F.determinant() < 0) continue;

                            output.push_back( F );
                            numberOfConfigurations++;
                            std::cout << numberOfConfigurations << std::endl;
                            if (numberOfConfigurations == maxConfigurations)
                                return output;
                        }
                    }
                }
            }
        }
        // sort the angles if rotations are being asked
        if (abs(maxStrain) < epsilon)
            std::transform(temp.begin(), temp.end(),
                           std::back_inserter(output),
                           [](const std::pair<IntScalarType,MatrixDimD>& p) {
                               return p.second;
                           });
        return output;
    }

    template<int dim> template<int dm>
    typename std::enable_if<dm==3,std::vector<LatticeVector<dim>>>::type
    Lattice<dim>::box(const std::vector<LatticeVector<dim>>& boxVectors, const std::string& filename) const
    {
        for(const LatticeVector<dim>& boxVector : boxVectors)
        {
            assert(this == &boxVector.lattice && "Box vectors belong to different lattice.");
        }

        // Form the box lattice
        MatrixDimD C;
        for (int i=0; i<dim; ++i) {
            C.col(i) = boxVectors[i].cartesian();
        }
        assert(C.determinant() != 0);
        Lattice<dim> boxLattice(C);

        auto planeParallelBasis= planeParallelLatticeBasis(boxVectors[1].cross(boxVectors[2]),true);
        MatrixDimD A;
        MatrixDimI mat;
        for(int i=0; i<dim; ++i) {
            A.col(i) = planeParallelBasis[i].cartesian();
            mat.col(i)=planeParallelBasis[i].latticeVector();
        }

        // l - lattice with plane parallel basis, where the second and third basis vectors
        //     lie on the plane spanned by boxVectors[1] and boxVectors[2]
        Lattice<dim> l(A);

        // Form plane parallel vectors v1 and v2
        // v1 is along boxVector[1]
        LatticeDirection<dim> v1(MatrixDimIExt<IntScalarType,dim>::adjoint(mat)*boxVectors[1],l);
        assert(v1.latticeVector()(0)==0);
        IntScalarType x,y;
        IntegerMath<IntScalarType>::extended_gcd(v1.latticeVector()(1),-v1.latticeVector()(2),x,y);
        LatticeVector<dim> temp(l);
        temp << 0,y,x;
        LatticeDirection<dim> v2(temp);

        int areaRatio= round((boxLattice.latticeBasis.col(1).cross(boxLattice.latticeBasis.col(2))).norm()/
                             (l.latticeBasis.col(1).cross(l.latticeBasis.col(2))).norm());
        int scale1= abs(IntegerMath<IntScalarType>::gcd(boxVectors[1]));
        int scale2= round(areaRatio/scale1);
        int scale0= round(abs(C.determinant()/latticeBasis.determinant()/areaRatio));

        std::vector<LatticeVector<dim>> output;

        // r0, r1, and r2 are reciprocal vectors perpendicular to areas spanned by boxVectors[0],
        // boxVectors[1], and boxVectors[2]
        ReciprocalLatticeVector<dim> r0_temp(*this), r1_temp(*this), r2_temp(*this);
        r0_temp= (boxVectors[1].cross(boxVectors[2])).reciprocalLatticeVector();
        if (r0_temp.dot(boxVectors[0]) < 0)
            r0_temp=-1*r0_temp;
        r1_temp= (boxVectors[2].cross(boxVectors[0])).reciprocalLatticeVector();
        if (r1_temp.dot(boxVectors[1]) < 0)
            r1_temp=-1*r1_temp;
        r2_temp= (boxVectors[0].cross(boxVectors[1])).reciprocalLatticeVector();
        if (r2_temp.dot(boxVectors[2]) < 0)
            r2_temp=-1*r2_temp;
        assert(r0_temp.dot(boxVectors[1])==0 && r0_temp.dot(boxVectors[2])==0);
        assert(r1_temp.dot(boxVectors[0])==0 && r1_temp.dot(boxVectors[2])==0);
        assert(r2_temp.dot(boxVectors[0])==0 && r2_temp.dot(boxVectors[1])==0);
        ReciprocalLatticeDirection<dim> r0(r0_temp), r1(r1_temp), r2(r2_temp);

        for(int i=0; i<scale0; ++i)
        {
            for(int j=0; j<scale1; ++j)
            {
                for(int k=0; k<scale2; ++k) {
                    LatticeVector<dim> vectorIn_l(l);
                    vectorIn_l << i, 0, 0;
                    vectorIn_l= vectorIn_l + j*v1.latticeVector()+k*v2.latticeVector();
                    //LatticeVector<dim> vector(latticeVector(vectorIn_l.cartesian()));
                    LatticeVector<dim> vector((mat*vectorIn_l).eval(),*this);
                    //LatticeDirection<dim> v1(MatrixDimIExt<IntScalarType,dim>::adjoint(mat)*boxVectors[1],l);

                    int c0= IntegerMath<IntScalarType>::positive_modulo(vector.dot(r0),boxVectors[0].dot(r0)) -
                            vector.dot(r0);
                    c0= c0/boxVectors[0].dot(r0);
                    int c1= IntegerMath<IntScalarType>::positive_modulo(vector.dot(r1),boxVectors[1].dot(r1)) -
                            vector.dot(r1);
                    c1= c1/boxVectors[1].dot(r1);
                    int c2= IntegerMath<IntScalarType>::positive_modulo(vector.dot(r2),boxVectors[2].dot(r2)) -
                            vector.dot(r2);
                    c2= c2/boxVectors[2].dot(r2);
                    vector= vector+c0*boxVectors[0]+c1*boxVectors[1]+c2*boxVectors[2];
                    output.push_back(vector);
                }
            }
        }

        if(!filename.empty()) {
            std::ofstream file;
            file.open(filename);
            if (!file) std::cerr << "Unable to open file";
            file << output.size() << std::endl;
            file << "Lattice=\"";

            for(const auto& vector:boxVectors) {
                file << vector.cartesian().transpose() << " ";
            }
            file << "\" Properties=atom_types:I:1:pos:R:3" << std::endl;

            for (const auto &vector: output)
                file << 1 << " " << vector.cartesian().transpose() << std::endl;

            file.close();
        }
        return output;
    }


    template<int dim> template<int dm>
    typename std::enable_if<dm==2,std::vector<LatticeVector<dim>>>::type
    Lattice<dim>::box(const std::vector<LatticeVector<dim>>& boxVectors, const std::string& filename) const
    {
        for(const LatticeVector<dim>& boxVector : boxVectors)
        {
            assert(this == &boxVector.lattice && "Box vectors belong to different lattice.");
        }

        // Form the box lattice
        MatrixDimD C;
        for (int i=0; i<dim; ++i) {
            C.col(i) = boxVectors[i].cartesian();
        }
        assert(C.determinant() != 0);
        Lattice<dim> boxLattice(C);

        // Form plane parallel vectors v0 and v1
        //LatticeDirection<dim> v0(this->latticeDirection(boxVectors[0].cartesian()));
        LatticeDirection<dim> v0(
                (VectorDimI) boxVectors[0]/IntegerMath<IntScalarType>::gcd(boxVectors[0]),
                *this);

        ReciprocalLatticeVector<dim> temp(*this);
        temp << v0.latticeVector()(1),-v0.latticeVector()(0);
        LatticeDirection<dim> v1(planeParallelLatticeBasis(temp,true)[0]);

        /*
        IntScalarType x,y;
        IntegerMath<IntScalarType>::extended_gcd(v0.latticeVector()(0),-v0.latticeVector()(1),x,y);
        LatticeVector<dim> temp(*this);
        temp << y,x;

        LatticeDirection<dim> v1(temp);
         */

        int areaRatio= round(abs(boxLattice.latticeBasis.determinant()/
                                 (*this).latticeBasis.determinant()));

        int scale1= abs(IntegerMath<IntScalarType>::gcd(boxVectors[0]));
        int scale2= round(areaRatio/scale1);

        std::vector<LatticeVector<dim>> output;

        // r0 and r1 are reciprocal vectors perpendicular to boxVectors[1] and boxVectors[0], respectively.
        ReciprocalLatticeVector<dim> r1_temp(*this), r0_temp(*this);
        r1_temp << boxVectors[0](1),-boxVectors[0](0);
        if (r1_temp.dot(boxVectors[1]) < 0)
            r1_temp=-1*r1_temp;
        r0_temp << boxVectors[1](1),-boxVectors[1](0);
        if (r0_temp.dot(boxVectors[0]) < 0)
            r0_temp=-1*r0_temp;
        assert(r1_temp.dot(boxVectors[0])==0);
        assert(r0_temp.dot(boxVectors[1])==0);
        ReciprocalLatticeDirection<dim> r1(r1_temp), r0(r0_temp);

        for(int j=0; j<scale1; ++j)
        {
            for(int k=0; k<scale2; ++k) {
                LatticeVector<dim> vector(j*v0.latticeVector()+k*v1.latticeVector());
                int c1= IntegerMath<IntScalarType>::positive_modulo(vector.dot(r1),boxVectors[1].dot(r1)) -
                        vector.dot(r1);
                c1= c1/boxVectors[1].dot(r1);
                int c2= IntegerMath<IntScalarType>::positive_modulo(vector.dot(r0),boxVectors[0].dot(r0)) -
                        vector.dot(r0);
                c2= c2/boxVectors[0].dot(r0);
                vector= vector+c1*boxVectors[1]+c2*boxVectors[0];
                output.push_back(vector);
            }
        }

        if(!filename.empty()) {
            std::ofstream file;
            file.open(filename);
            if (!file) std::cerr << "Unable to open file";
            file << output.size() << std::endl;
            file << "Lattice=\"";

            for(const auto& vector:boxVectors) {
                file << vector.cartesian().transpose() << " 0 ";
            }
            file << " 0 0 1 ";
            file << "\" Properties=atom_types:I:1:pos:R:3" << std::endl;

            for (const auto &vector: output)
                file << 1 << " " << vector.cartesian().transpose() << " 0 " << std::endl;

            file.close();
        }
        return output;
    }


    template class Lattice<1>;

    template class Lattice<2>;
    template std::vector<typename Lattice<2>::MatrixDimD> Lattice<2>::generateCoincidentLattices<2>(
            const double &maxStrain, const int &maxDen, const int &N) const;
    template std::vector<typename Lattice<2>::MatrixDimD> Lattice<2>::generateCoincidentLattices<2>(
            const Lattice<2> &undeformedLattice, const double &maxStrain, const int &maxDen, const int &N) const;
    template std::vector<LatticeVector<2>> Lattice<2>::box<2>(const std::vector<LatticeVector<2>> &boxVectors,
                                                           const std::string &filename) const;


    template class Lattice<3>;
    template std::vector<typename Lattice<3>::MatrixDimD> Lattice<3>::generateCoincidentLattices<3>(
            const ReciprocalLatticeDirection<3> &rd, const double &maxDen, const int& N) const;
    template std::vector<LatticeVector<3>> Lattice<3>::box<3>(const std::vector<LatticeVector<3>> &boxVectors,
                                                              const std::string &filename) const;

    template class Lattice<4>;
    template class Lattice<5>;
}
#endif
