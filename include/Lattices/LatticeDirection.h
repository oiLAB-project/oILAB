/* This file is part of gbLAB.
 *
 * gbLAB is distributed without any warranty under the MIT License.
 */


#ifndef gbLAB_LatticeDirection_h_
#define gbLAB_LatticeDirection_h_

#include <LatticeModule.h>

namespace gbLAB
{
    template <int dim>
    struct LatticeDirection : public LatticeVector<dim>
    {
        typedef typename LatticeCore<dim>::IntScalarType IntScalarType;
        
        LatticeDirection(const LatticeVector<dim>& v) ;

        LatticeDirection(const Eigen::Matrix<IntScalarType,dim,1>& v,
                         const Lattice<dim>& lat) ;
        
    };
    
} // end namespace
#endif
