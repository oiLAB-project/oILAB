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
    /* inherits */ protected ReciprocalLatticeVector<dim>
    {
        typedef typename LatticeCore<dim>::IntScalarType IntScalarType;


        ReciprocalLatticeDirection(const ReciprocalLatticeDirection<dim>& other) = default;
        ReciprocalLatticeDirection(const ReciprocalLatticeVector<dim>& v) ;

        using ReciprocalLatticeVector<dim>::cartesian;
        using ReciprocalLatticeVector<dim>::lattice;
        using ReciprocalLatticeVector<dim>::dot;

        const ReciprocalLatticeVector<dim>& reciprocalLatticeVector() const
        {
            return static_cast<const ReciprocalLatticeVector<dim>&>(*this);
        }

    };

    template<int dim>
    basic_ostream<char>& operator<<(basic_ostream<char>& s, const ReciprocalLatticeDirection<dim>& m);
} // end namespace
#endif
