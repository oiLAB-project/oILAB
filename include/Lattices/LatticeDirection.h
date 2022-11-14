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
    struct LatticeDirection : protected LatticeVector<dim>
    {
        using IntScalarType = typename LatticeCore<dim>::IntScalarType ;
        using VectorDimI = typename LatticeCore<dim>::VectorDimI;

        
        LatticeDirection(const LatticeVector<dim>& v);
        LatticeDirection(const LatticeDirection<dim>& other) = default;

        LatticeDirection(const VectorDimI& v,
                         const Lattice<dim>& lat) ;

        using LatticeVector<dim>::cartesian;
        using LatticeVector<dim>::lattice;
        using LatticeVector<dim>::dot;

        const LatticeVector<dim>& latticeVector() const
        {
            return static_cast<const LatticeVector<dim>&>(*this);
        }
    };

    template<int dim>
    basic_ostream<char> &operator<<(basic_ostream<char> &s, const LatticeDirection<dim>& m);


} // end namespace
#endif
