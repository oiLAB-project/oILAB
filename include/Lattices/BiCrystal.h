/* This file is part of gbLAB.
 *
 * gbLAB is distributed without any warranty under the MIT License.
 */


#ifndef gbLAB_BiCrystal_h_
#define gbLAB_BiCrystal_h_

#include <LatticeModule.h>
#include <SmithDecomposition.h>
#include "RationalMatrix.h"
#include "LLL.h"
#include "RLLL.h"
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
         * Outputs lattice vector in lattice \f$\mathcal D\f$ that is equal to the inputted vector \f$\textbf v\f$
         * that belongs to \f$\mathcal A\f$, \f$\mathcal B\f$, \f$\mathcal C\f$, or \f$\mathcal D\f$
         * @param v - lattice vector
         * @return LatticeVector in \f$\mathcal D\f$
         */
        LatticeVector<dim> getLatticeVectorInD(const LatticeVector<dim>& v) const;
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
        typename std::enable_if<dm==2 || dm==3,std::map<IntScalarType,Gb<dm>>>::type
        generateGrainBoundaries(const LatticeDirection<dim>& d, int div=30) const;


        /*! This function outputs/prints a 2D bicrystal (two lattices that form the GB and
         * the CSL) bounded by a box defined using
         * two input box vectors. The box vectors have to be linearly independent lattice
         * vectors. The function optimizes boxVectors[0]
         * to make the box as orthogonal as possible depending on the \p orthogonality parameter.
         *
         *
         * @tparam dm dimension (int)
         * @param boxVectors two linearly independent lattice vectors.
         * @param orthogonality (double) a value in the interval \f$[0,1]\f$.
         * @param filename (optional) name of the output file
         * @param orient (optional) While printing to a file, orient the system such that one of the box sides
         * is along the global x axis. This flag does not
         * influence the returning configuration, only the configuration printed to the file.
         * @return lattice points of the bicrystal (along with the CSL) bounded by the box (std::vector<LatticeVector<2>>).
         */
        template<int dm=dim>
        typename std::enable_if<dm==2 || dm==3,std::vector<LatticeVector<dim>>>::type
        box(std::vector<LatticeVector<dim>>& boxVectors, 
                const double& orthogonality, 
                const int& dsclFactor,
                std::string filename= "", 
                bool orient=false) const;
    };
    
    
} // end namespace
#endif

