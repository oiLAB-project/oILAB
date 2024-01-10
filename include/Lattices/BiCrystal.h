/* This file is part of gbLAB.
 *
 * gbLAB is distributed without any warranty under the MIT License.
 */


#ifndef gbLAB_BiCrystal_h_
#define gbLAB_BiCrystal_h_

#include <LatticeModule.h>
#include <SmithDecomposition.h>
#include <RationalMatrix.h>
#include <LLL.h>
#include <RLLL.h>
#include <unordered_set>
#include "Rotation.h"


namespace gbLAB
{
    /*!Class template that computes the coincident-site-lattice (CSL) of two
     * parent lattices using the Smith Normal Form [1].
     *
     * [1] Coincidence Lattices and Associated Shear Transformations
     */
    template <int dim>
    class BiCrystal : public RationalMatrix<dim>
    /*             */,public SmithDecomposition<dim>
    {
        typedef typename LatticeCore<dim>::IntScalarType IntScalarType;
        typedef typename LatticeCore<dim>::VectorDimD VectorDimD;
        typedef typename LatticeCore<dim>::MatrixDimD MatrixDimD;
        typedef typename LatticeCore<dim>::VectorDimI VectorDimI;
        typedef typename LatticeCore<dim>::MatrixDimI MatrixDimI;

        
        static MatrixDimI getM(const RationalMatrix<dim>& rm, const SmithDecomposition<dim>& sd);
        static MatrixDimI getN(const RationalMatrix<dim>& rm, const SmithDecomposition<dim>& sd);
        static MatrixDimI getLambdaA(const MatrixDimI& M, const MatrixDimI& N);
        static MatrixDimI getLambdaB(const MatrixDimI& M, const MatrixDimI& N);
        static MatrixDimD getCSLBasis(const Lattice<dim>& A,
                                      const Lattice<dim>& B,
                                      const SmithDecomposition<dim>& sd,
                                      const MatrixDimI& M,
                                      const MatrixDimI& N,
                                      const bool& useRLLL);
        static MatrixDimD getDSCLBasis(const Lattice<dim>& A,
                                       const Lattice<dim>& B,
                                       const SmithDecomposition<dim>& sd,
                                       const MatrixDimI& M,
                                       const MatrixDimI& N,
                                       const bool& useRLLL);

    public:

        const Lattice<dim>& A;
        const Lattice<dim>& B;

        /*! \brief Integer matrix that connects the bases \f$ \textbf A^{\|} \f$ and \f$ \textbf C^{\|} \f$ of
         * lattices \f$\mathcal A\f$ and the CSL \f$\mathcal C\f$, respectively: \f$ \textbf C^{\|} = \textbf A^{\|} \textbf M\f$.
         */
        const MatrixDimI M;

        /*! \brief Integer matrix that connects the bases \f$ \textbf B^{\|} \f$ and \f$ \textbf C^{\|} \f$ of
         * lattices \f$\mathcal B\f$ and the CSL \f$\mathcal C\f$, respectively: \f$ \textbf C^{\|} = \textbf B^{\|} \textbf N\f$.
         */
        const MatrixDimI N;

        /*! \brief Signed ratio of the unit cell volume of \f$\mathcal C\f$ to that of \f$\mathcal A\f$.
         *  \f$ \Sigma_{\mathcal A} = \det(\textbf M)\f$.
         */
        const IntScalarType sigmaA;

        /*! \brief Signed ratio of the unit cell volume of \f$\mathcal C\f$ to that of \f$\mathcal B\f$:
         *  \f$ \Sigma_{\mathcal B} = \det(\textbf N)\f$.
         */
        const IntScalarType sigmaB;

        /*! \brief \f$ \Sigma = |\Sigma_{\mathcal A}|\f$ if \f$ \Sigma_{\mathcal A} = \Sigma_{\mathcal B}\f$, else
         *  \f$ \Sigma = 0 \f$
         */
        const int sigma;

        /*! \brief CSL lattice \f$\mathcal C\f$.
         */
        const Lattice<dim> csl;

        /*! \brief DCSL lattice \f$\mathcal D\f$.
         */
        const Lattice<dim> dscl;

        /*! \brief Lattice \f$\mathcal A\f$ with basis \f$\textbf A^\|\f$
         */
        const Lattice<dim> Ap;

        /*! \brief Lattice \f$\mathcal B\f$ with basis \f$\textbf B^\|\f$
         */
        const Lattice<dim> Bp;

        /*! \brief Shift tensor \f$\Lambda_{\mathcal A}:\mathbb Z_{\mathcal D} \to \mathbb Z_{\mathcal D}\f$ describes the
         * shift in the CSL when lattice \f$\mathcal A\f$ is shifted by a DSCL vector.
         */
        const MatrixDimI LambdaA;

        /*! \brief Shift tensor \f$\Lambda_{\mathcal B}:\mathbb Z_{\mathcal D} \to \mathbb Z_{\mathcal D}\f$ describes the
         * shift in the CSL when lattice \f$\mathcal B\f$ is shifted by a DSCL vector.
         */
        const MatrixDimI LambdaB;


        LatticeVector<dim> shiftTensorA(const LatticeVector<dim>& d) const;
        LatticeVector<dim> shiftTensorB(const LatticeVector<dim>& d) const;

        /**********************************************************************/
        /*! \brief Constructs a bicrystal from two lattices \f$\mathcal A \f$ and \f$\mathcal B \f$ by computing the
         *  parallel bases \f$\textbf A^\|, \textbf B^\|, \textbf C^\|, \textbf D^\|\f$ for lattices
         *  \f$\mathcal A, \mathcal B, \mathcal C, \mathcal D\f$, respectively.
         *  If the flag useRLLL is .true., then the bases of CSL and DSCL are reduced using the LLL algorithm.
         *
         * \param[in] Lattices  \f$\mathcal A \f$ and \f$\mathcal B \f$, and useRLLL flag
         * \returns   A bicrystal object
         * */
        BiCrystal(const Lattice<dim>& A,
                  const Lattice<dim>& B,
                  const bool& useRLLL=false);

        /*!
         * Outputs lattice vector in lattice \f$\mathcal A\f$ that is equal to the inputted vector \f$\textbf v\f$
         * that belongs to \f$\mathcal A\f$ or \f$\mathcal C\f$
         * @param v - lattice vector
         * @return LatticeVector in \f$\mathcal A\f$
         */
        LatticeVector<dim> getLatticeVectorInA(const LatticeVector<dim>& v) const;
        /*!
         * Outputs lattice vector in lattice \f$\mathcal B\f$ that is equal to the inputted vector \f$\textbf v\f$
         * that belongs to \f$\mathcal B\f$ or \f$\mathcal C\f$
         * @param v - lattice vector
         * @return LatticeVector in \f$\mathcal B\f$
         */
        LatticeVector<dim> getLatticeVectorInB(const LatticeVector<dim>& v) const;
        /*!
         * Outputs lattice direction in the CSL \f$\mathcal C\f$ that is parallel to the inputted vector \f$\textbf v\f$
         * that belongs to one of the four lattices, \f$\mathcal A\f$, \f$\mathcal B\f$, \f$\mathcal C\f$, or \f$\mathcal D\f$,
         * @param v - lattice vector
         * @return LatticeDirection in \f$\mathcal C\f$
         */
        LatticeDirection<dim> getLatticeDirectionInC(const LatticeVector<dim>& v) const;
        /*!
         * Outputs lattice direction in the DSCL \f$\mathcal D\f$ that is parallel to the inputted vector \f$\textbf v\f$
         * that belongs to one of the four lattices, \f$\mathcal A\f$, \f$\mathcal B\f$, \f$\mathcal C\f$, or \f$\mathcal D\f$,
         * @param v - lattice vector
         * @return LatticeDirection in \f$\mathcal D\f$
         */
        LatticeDirection<dim> getLatticeDirectionInD(const LatticeVector<dim>& v) const;

        /*!
         * Outputs reciprocal lattice direction in the dual lattice \f$\mathcal A^*\f$ that is parallel to the inputted
         * reciprocal vector \f$\textbf v\f$ that belongs to one of the four dual lattices, \f$\mathcal A^*\f$, \f$\mathcal B^*\f$,
         * \f$\mathcal C^*\f$, or \f$\mathcal D^*\f$,
         * @param v - lattice vector
         * @return ReciprocalLatticeDirection in \f$\mathcal A^*\f$
         */
        ReciprocalLatticeDirection<dim> getReciprocalLatticeDirectionInA(const ReciprocalLatticeVector<dim>& v) const;
        /*!
         * Outputs reciprocal lattice direction in the dual lattice \f$\mathcal B^*\f$ that is parallel to the inputted
         * reciprocal vector \f$\textbf v\f$ that belongs to one of the four dual lattices, \f$\mathcal A^*\f$, \f$\mathcal B^*\f$,
         * \f$\mathcal C^*\f$, or \f$\mathcal D^*\f$,
         * @param v - lattice vector
         * @return ReciprocalLatticeDirection in \f$\mathcal B^*\f$
         */
        ReciprocalLatticeDirection<dim> getReciprocalLatticeDirectionInB(const ReciprocalLatticeVector<dim>& v) const;
        /*!
         * Outputs reciprocal lattice direction in the dual lattice \f$\mathcal C^*\f$ that is parallel to the inputted
         * reciprocal vector \f$\textbf v\f$ that belongs to one of the four dual lattices, \f$\mathcal A^*\f$, \f$\mathcal B^*\f$,
         * \f$\mathcal C^*\f$, or \f$\mathcal D^*\f$,
         * @param v - lattice vector
         * @return ReciprocalLatticeDirection in \f$\mathcal C^*\f$
         */
        ReciprocalLatticeDirection<dim> getReciprocalLatticeDirectionInC(const ReciprocalLatticeVector<dim>& v) const;
        /*!
         * Outputs reciprocal lattice direction in the dual lattice \f$\mathcal D^*\f$ that is parallel to the inputted
         * reciprocal vector \f$\textbf v\f$ that belongs to one of the four dual lattices, \f$\mathcal A^*\f$, \f$\mathcal B^*\f$,
         * \f$\mathcal C^*\f$, or \f$\mathcal D^*\f$,
         * @param v - lattice vector
         * @return ReciprocalLatticeDirection in \f$\mathcal D^*\f$
         */
        ReciprocalLatticeDirection<dim> getReciprocalLatticeDirectionInD(const ReciprocalLatticeVector<dim>& v) const;

        /*!
         * \brief Given a tilt axis \f$\textbf d\f$, that belongs to lattices \f$\mathcal A\f$ or \f$\mathcal B\f$, this
         * function generate a set of tilt GBs. CURRENTLY ONLY WORDS FOR DIMENSION 3
         * @tparam dm
         * @param d - LatticeDirection that describes the tilt axis
         * @param div - parameter to span the GBs
         * @return A data structure that stores GBs sorted in increasing order of their inclination angle.
         */
        template<int dm=dim>
        typename std::enable_if<dm==2 || dim==3,std::map<IntScalarType,Gb<dm>>>::type
        generateGrainBoundaries(const LatticeDirection<dim>& d, int div=30) const
        {
            if (&d.lattice != &A && &d.lattice != &B)
                throw std::runtime_error("The tilt axis does not belong to lattices A and B  \n");
            std::vector<Gb<dm>> gbVec;
            double epsilon=1e-8;
            bool stgbExists= false;
            int count= -1;
            int stgbCount= 0;
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
                if (stgbExists)
                    cosAngle= gb.nA.cartesian().normalized().dot(gbVec[stgbCount].nA.cartesian().normalized());
                else
                    cosAngle= gb.nA.cartesian().normalized().dot(gbVec[0].nA.cartesian().normalized());
                if (cosAngle-1>-epsilon) cosAngle= 1.0;
                if (cosAngle+1<epsilon) cosAngle= -1.0;

                double angle= acos(cosAngle);
                IntScalarType key= angle*keyScale;
                gbSet.insert(std::pair<IntScalarType,Gb<dm>>(key,gb));
            }
            return gbSet;
        }


        template<int dm=dim>
        typename std::enable_if<dm==2,std::vector<LatticeVector<dim>>>::type
        box(std::vector<LatticeVector<dim>> boxVectors, double const& orthogonality, std::string filename= "", bool orient=false) const
        {
            assert(orthogonality>=0.0 && orthogonality<=1.0 &&
                           "The \"orthogonality\" parameter should be between 0.0 and 1.0");
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
            assert(C.determinant() != 0);


            // Adjust boxVector[0] such that it is as orthogonal as possible to boxVector[1]
            auto boxVectorTemp= boxVectors[0];
            ReciprocalLatticeVector<dim> temp(csl);
            temp << -boxVectors[1](1),boxVectors[1](0);
            ReciprocalLatticeDirection<dim> nC(temp);
            std::cout << "Input box height= " << abs(nC.planeSpacing()*boxVectors[0].dot(nC)) << std::endl;
            auto basis= csl.planeParallelLatticeBasis(nC,true);

            int planesToExplore= nC.stacking();
            std::cout << "Exploring " << planesToExplore << " planes" << std::endl;
            MatrixDimI boxLatticeIndices;
            boxLatticeIndices.col(0)= boxVectors[0];
            for (int i=1; i<dim; ++i)
                boxLatticeIndices.col(i)= boxVectors[i]/IntegerMath<long long int>::gcd(boxVectors[i]);
            double minAngle= M_PI/2;
            double nonOrthogonality= 1;
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
                double angle= abs(acos(boxVectorTemp.cartesian().normalized().dot(rC.cartesian().normalized())));
                if(angle < minAngle) {
                    minAngle= angle;
                    boxVectorUpdated= boxVectorTemp;
                    if (angle < (1-orthogonality)*M_PI/2)
                        break;
                }

            }
            boxVectors[0]=boxVectorUpdated;
            C.col(0)= boxVectors[0].cartesian();
            std::cout << "Updated box height= " << abs(nC.planeSpacing()*boxVectors[0].dot(nC)) << std::endl;


            MatrixDimD rotation= Eigen::Matrix<double,dim,dim>::Identity();;
            if (orient) {
                Eigen::Matrix<double,dim,dim-1> orthogonalVectors;
                orthogonalVectors.col(0)= C.col(1).normalized();
                rotation = Rotation<dim>(orthogonalVectors);
            }
            assert((rotation*rotation.transpose()).template isApprox(Eigen::Matrix<double,dim,dim>::Identity())
                   && "Cannot orient the grain boundary. Box vectors are not orthogonal.");

            std::vector<LatticeVector<dim>> configurationA, configurationB, configurationC, configurationD;
            std::vector<LatticeVector<dim>> configuration;

            std::vector<LatticeVector<dim>> boxVectorsInA, boxVectorsInB, boxVectorsInD;
            for(const auto& boxVector : boxVectors) {
                boxVectorsInA.push_back(getLatticeVectorInA(boxVector));
                boxVectorsInB.push_back(getLatticeVectorInB(boxVector));
            }

            configurationA= A.box(boxVectorsInA);
            configurationB= B.box(boxVectorsInB);
            configurationC= csl.box(boxVectors);

            configuration= configurationA;
            configuration.insert(configuration.end(),configurationB.begin(),configurationB.end());
            std::cout << "Number of lattice points in A = " << configurationA.size() << std::endl;
            std::cout << "Number of lattice points in B = " << configurationB.size() << std::endl;
            std::cout << "Number of lattice points in CSL = " << configurationC.size() << std::endl;

            if(!filename.empty()) {
                std::ofstream file;
                file.open(filename);
                if (!file) std::cerr << "Unable to open file";
                file << configurationA.size() + configurationB.size() + configurationC.size() << std::endl;
                file << "Lattice=\"";

                if (dim==2) {
                    file << (rotation*boxVectors[0].cartesian()).transpose() << " 0 ";
                    file << (rotation*boxVectors[1].cartesian()).transpose() << " 0 ";
                    file << " 0 0 1 ";
                    file << "\" Properties=atom_types:I:1:pos:R:3:radius:R:1" << std::endl;
                    for (const auto &vector: configurationA)
                        file << 1 << " " << (rotation*vector.cartesian()).transpose() << " " << 0.0 << "  " << 0.05 << std::endl;
                    for (const auto &vector: configurationB)
                        file << 2 << " " << (rotation*vector.cartesian()).transpose() << " " << 0.0 << "  " << 0.05 << std::endl;
                    for (const auto &vector: configurationC)
                        file << 3 << " " << (rotation*vector.cartesian()).transpose() << " " << 0.0 << "  " << 0.2 << std::endl;
                }
                else if (dim==3){
                    file << (rotation*boxVectors[0].cartesian()).transpose()  << " ";
                    file << (rotation*boxVectors[1].cartesian()).transpose() << " ";
                    file << (rotation*boxVectors[2].cartesian()).transpose() << " ";
                    file << "\" Properties=atom_types:I:1:pos:R:3:radius:R:1" << std::endl;

                    for (const auto &vector: configurationA)
                        file << 1 << " " << (rotation*vector.cartesian()).transpose() << "  " << 0.05 << std::endl;
                    for (const auto &vector: configurationB)
                        file << 2 << " " << (rotation*vector.cartesian()).transpose() << "  " << 0.05 << std::endl;
                    for (const auto &vector: configurationC)
                        file << 3 << " " << (rotation*vector.cartesian()).transpose() << "  " << 0.2 << std::endl;
                }


                file.close();
            }
            return configuration;
        }
    };
    
    
} // end namespace
#endif

