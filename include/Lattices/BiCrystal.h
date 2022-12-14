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
        const int sigmaA;

        /*! \brief Signed ratio of the unit cell volume of \f$\mathcal C\f$ to that of \f$\mathcal B\f$:
         *  \f$ \Sigma_{\mathcal B} = \det(\textbf N)\f$.
         */
        const int sigmaB;

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

        LatticeDirection<dim> getLatticeDirectionInC(const LatticeVector<dim>& v) const;
        LatticeDirection<dim> getLatticeDirectionInD(const LatticeVector<dim>& v) const;

        ReciprocalLatticeDirection<dim> getReciprocalLatticeDirectionInA(const ReciprocalLatticeVector<dim>& v) const;
        ReciprocalLatticeDirection<dim> getReciprocalLatticeDirectionInB(const ReciprocalLatticeVector<dim>& v) const;
        ReciprocalLatticeDirection<dim> getReciprocalLatticeDirectionInC(const ReciprocalLatticeVector<dim>& v) const;
        ReciprocalLatticeDirection<dim> getReciprocalLatticeDirectionInD(const ReciprocalLatticeVector<dim>& v) const;

        template<int dm=dim>
        typename std::enable_if<dm==2 || dim==3,std::map<IntScalarType,Gb<dm>>>::type
        generateGrainBoundaries(const LatticeDirection<dim>& d, int div=30) const
        {
            // asert that d should be a lattice vector of A or B
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

                            VectorDimI nAintegerCoords = gb.nA.reciprocalLatticeVector();
                            VectorDimI nBintegerCoords = gb.nB.reciprocalLatticeVector();

                            nAintegerCoords = nAintegerCoords.cwiseAbs();
                            nBintegerCoords = nBintegerCoords.cwiseAbs();
                            std::sort(nAintegerCoords.data(), nAintegerCoords.data() + dm);
                            std::sort(nBintegerCoords.data(), nBintegerCoords.data() + dm);
                            if (nAintegerCoords == nBintegerCoords) {
                                stgbExists = true;
                                stgbCount = count;
                            }
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
    };
    
    
} // end namespace
#endif

