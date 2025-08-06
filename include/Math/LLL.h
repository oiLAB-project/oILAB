/* This file is part of gbLAB.
 *
 * gbLAB is distributed without any warranty under the MIT License.
 */


#ifndef gbLAB_LLL_h_
#define gbLAB_LLL_h_

#include <cmath>
#include <algorithm> // std::max
#include <utility>   // std::swap
#include "Eigen/Dense"
#include <iostream>
// http://www.arageli.org/download
// https://www.mathworks.com/matlabcentral/fileexchange/49457-lattice-reduction-mimo?focused=3859922&tab=function


namespace gbLAB
{
    /*! An implementation of the  Lenstra–Lenstra–Lovász (LLL) lattice basis
     *  reduction algorithm for integers.
     */
    class LLL
    {
        
        
        
        /**********************************************************************/
        void lll_gram_schmidt_int(const int& k);
        
        /**********************************************************************/
        void lll_size_reduction_int(const int& k,
                                    const int& l);
        
        /**********************************************************************/
        void lll_interchange_int(const int&  k,
                                 const int&  k_max);
        
        Eigen::MatrixXi B;
        Eigen::VectorXi  d;
        Eigen::MatrixXi  H;
        Eigen::MatrixXi  Lambda;
        
        
    public:
        
        template<int m,int n>
        LLL(const Eigen::Matrix<int,m,n>& B_in) ;
    };
    
} // end namespace
#endif
