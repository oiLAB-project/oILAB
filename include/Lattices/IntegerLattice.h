/* This file is part of gbLAB.
 *
 * gbLAB is distributed without any warranty under the MIT License.
 */


#ifndef gbLAB_IntegerLattice_h_
#define gbLAB_IntegerLattice_h_

#include <DiophantineSolver.h>

namespace gbLAB
{
    /*! \brief Lattice class
     *
     *  lattice class description
     * */
    template <int dim>
    class IntegerLattice
    {
        using IntScalarType = long long int;
        using VectorDimI = Eigen::Matrix<IntScalarType,dim,1>;
        using MatrixDimI = Eigen::Matrix<IntScalarType,dim,dim>;
    public:
        const MatrixDimI    latticeBasis;
        IntegerLattice() : latticeBasis(MatrixDimI::Identity())
        {}

        static Eigen::Matrix<IntScalarType,dim-1,dim>
        perpendicularDirections(const Eigen::Vector<IntScalarType,dim>& in)
        {
            Eigen::Matrix<IntScalarType,dim-1,dim> out= Eigen::Matrix<IntScalarType,dim-1,dim>::Zero();

            int firstNonZero= 0;
            for (int ind= 0; ind<dim; ind++)
            {
                if (in(ind) != 0)
                {
                    firstNonZero = ind;
                    break;
                }
                else
                {
                    out(ind,ind)= 1;
                }
            }

            if (firstNonZero==dim-1) return out;

            out(firstNonZero,firstNonZero)= -in(firstNonZero+1)/IntegerMath<IntScalarType>::gcd(in(firstNonZero),in(firstNonZero+1));
            out(firstNonZero,firstNonZero+1)= in(firstNonZero)/IntegerMath<IntScalarType>::gcd(in(firstNonZero),in(firstNonZero+1));

            for (int ind= firstNonZero+2; ind<dim; ind++)
            {
                if (in(ind)==0)
                {
                    out.row(ind-1)= Eigen::Vector<IntScalarType, dim>::Zero(dim);
                    out(ind-1,ind)= 1;
                    continue;
                }
                Eigen::Vector<IntScalarType,3> truncated_in;
                truncated_in << in(ind-2), in(ind-1), in(ind);
                //truncated_in= truncated_in/IntegerMath<IntScalarType>::gcd(truncated_in);

                IntScalarType x,y;
                DiophantineSolver<IntScalarType>::solveDiophantine2vars(truncated_in(0),
                                                         truncated_in(1),
                                                         -truncated_in(2)*IntegerMath<IntScalarType>::gcd(truncated_in(0),truncated_in(1)),x,y);

                out(ind-1,ind)= IntegerMath<IntScalarType>::gcd(truncated_in(0),truncated_in(1));
                out(ind-1,ind-2)= x;
                out(ind-1,ind-1)= y;

                out.row(ind-1)= out.row(ind-1)/IntegerMath<IntScalarType>::gcd(out.row(ind-1));
            }
            return out;
        }
    };
}
#endif
