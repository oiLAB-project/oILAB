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
         * lattices \f$\mathcal A\f$ and the CSL \f$\mathcal C\f$, respectively.
         */
        const MatrixDimI M;
        const MatrixDimI N;
        const int sigmaA;
        const int sigmaB;
        const int sigma;
        const Lattice<dim> csl;
        const Lattice<dim> dscl;
        const Lattice<dim> Ap;
        const Lattice<dim> Bp;
        
        /**********************************************************************/
        /*! \brief Constructs a bicrystal from two lattices \f$\textbf A \f$ and \f$\textbf B \f$ by computing the
         *  parallel bases Ap and Bp, CSL, and the DSCL.
         *  If the flass useRLLL is .true., then the bases of CSL and DSCL are reduced using the LLL algorithm.
         *
         * \param[in] Lattices  \f$\mathcal A \f$ and \f$\mathcal B \f$, and useRLLL flag
         * \returns   A bicrystal object
         * */
        BiCrystal(const Lattice<dim>& A,
                  const Lattice<dim>& B,
                  const bool& useRLLL=true);
        
//        LatticeDirection<dim> AtoCSLvector(const LatticeVector<dim>& v) const;
        
    };
    
    
} // end namespace
#endif

