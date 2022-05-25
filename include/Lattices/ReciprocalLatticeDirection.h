/* This file is part of gbLAB.
 *
 * gbLAB is distributed without any warranty under the MIT License.
 */


#ifndef gbLAB_ReciprocalLatticeDirection_h_
#define gbLAB_ReciprocalLatticeDirection_h_

#include <LatticeModule.h>

namespace gbLAB
{
    template <int dim>
    struct ReciprocalLatticeDirection :
    /* inherits */ public ReciprocalLatticeVector<dim>
    {
        typedef typename LatticeCore<dim>::IntScalarType IntScalarType;


        ReciprocalLatticeDirection(const ReciprocalLatticeVector<dim>& v) ;
                
    };
    
} // end namespace
#endif
