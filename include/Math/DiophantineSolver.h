/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef gbLAB_DiophantineSolver_h_
#define gbLAB_DiophantineSolver_h_

#include <Eigen/Dense>
#include <IntegerMath.h>
#include <climits>
#include <cmath>
#include <sortIndices.h>
namespace gbLAB
{
    template<typename IntScalarType>
    class DiophantineSolver
    {
        public:
        // Finds such x and y, that a  x + b  y = gcd(a, b)
        //https://github.com/ADJA/algos/blob/master/NumberTheory/DiophantineEquation.cpp
        static IntScalarType extended_gcd(IntScalarType a, IntScalarType b, IntScalarType &x, IntScalarType &y)
        {
            if (b == 0)
            {
                x = 1;
                y = 0;
                return a;
            }
            IntScalarType x1, y1;
            IntScalarType g = extended_gcd(b, a % b, x1, y1);
            x = y1;
            y = x1 - (a / b) * y1;
            return g;
        }
/*
        // Find u such that a_1 u_1 + a_2 u_2 + .... + a_n u_n = 1
        // method 1 - does not work, the algorithm is incorrect

        // "A fast algorithm to find reduced hyperplane unit cells and solve N-dimensional Bézout's identities."
        //  Acta Crystallographica Section A: Foundations and Advances 77.5 (2021).

        static Eigen::Vector<IntScalarType,Eigen::Dynamic> solveBezout(const Eigen::Vector<IntScalarType,Eigen::Dynamic>& p)
        {
            int n= p.size();
            Eigen::Vector<IntScalarType,Eigen::Dynamic> out(n);
            out= Eigen::Vector<IntScalarType,Eigen::Dynamic>::Zero(n);

            if (n==2)
            {
                IntScalarType g= extended_gcd(p(0),p(1),out(0),out(1));
                std::cout << p(0) << "   " << p(1) << std::endl;
                if (g<0) out= -out;
                return out;
            }

            int ind=0;
            for (auto pi: p)
            {
                if (abs(pi)==1)
                {
                    out(ind) = pi;
                    return out;
                }
                ind++;
            }

            Eigen::Vector<IntScalarType,Eigen::Dynamic> pS(n);
            ind= -1;
            for (auto i: sortIndices(p))
            {
                ind++;
                if (p(i) == 0) break;
                pS(ind)= p(i);
            }

            IntScalarType pSi0= pS(ind);
            int numNonZeros= ind+1;
            Eigen::Vector<IntScalarType,Eigen::Dynamic> quotients(numNonZeros-1);
            Eigen::Vector<IntScalarType,Eigen::Dynamic> remainders(numNonZeros-1);
            for(int i=0; i<quotients.size(); i++)
            {
                quotients(i)= pS(i)/pSi0;
                remainders(i)= pS(i)%pSi0;

                //if (pS(i)/pSi0 < 0)
                //    quotients(i)= floor((double)pS(i)/pSi0);
                //else
                //    quotients(i)= pS(i)/pSi0;
                //remainders(i)= pS(i)-pSi0*quotients(i);
            }

            Eigen::Vector<IntScalarType,Eigen::Dynamic> u(numNonZeros-1);
            u= solveBezout(remainders);

            out(Eigen::seq(0,numNonZeros-2))= u;
            out(numNonZeros-1)= -quotients.dot(u);

            Eigen::Vector<IntScalarType,Eigen::Dynamic> outRestored(n);
            ind= 0;
            for (auto i : sortIndices(p))
            {
                outRestored(i) = out(ind);
                ind++;
            }
            return outRestored;
        }
*/
        // Find u such that a_1 u_1 + a_2 u_2 + .... + a_n u_n = 1
        // Method 0:
        // "A fast algorithm to find reduced hyperplane unit cells and solve N-dimensional Bézout's identities."
        //  Acta Crystallographica Section A: Foundations and Advances 77.5 (2021).
        template<typename ArrayType>
        static ArrayType solveBezout(const ArrayType& a)
        {
            int n= a.size();
            if (n<2) throw std::runtime_error("the size of arrays should be at least two");
            if (abs(IntegerMath<IntScalarType>::gcd(a)) != 1)
                throw std::runtime_error("No solution since the gcd is not 1 or -1.");

            ArrayType u(n);

            IntScalarType p,q;
            IntScalarType g12= extended_gcd(a(0),a(1),p,q); // g12 - gcd of a(0) and a(1)

            // form the n-1 sized vector na= {g,a2,...,a_{n-1}}
            Eigen::Vector<IntScalarType,Eigen::Dynamic> na(n-1);
            na(0)= g12; na(Eigen::seq(1,n-2))= a(Eigen::seq(2,n-1));

            Eigen::Vector<IntScalarType,Eigen::Dynamic> k(n-1);
            if (n==2)
            {
                IntScalarType g= extended_gcd(a(0),a(1),u(0),u(1));
                if (g<0) u= -u;
                return u;
            }
            else
            {
                k = solveBezout(na);
                // {p*k0,q*k0,k1,..,k_{n-2}}
                u(0) = p * k(0);
                u(1) = q * k(0);
                u(Eigen::seq(2, n - 1)) = k(Eigen::seq(1, n - 2));
                return u;
            }
        }

        // Solves equation a  x + b  y = c, writes answer to x and y
        static void solveDiophantine2vars(IntScalarType a, IntScalarType b, IntScalarType c, IntScalarType &x, IntScalarType &y)
        {

            IntScalarType g = extended_gcd(a, b, x, y);

            if (c % g != 0)
            {
                std::cout << a << "  " << b << "  " << c << "   " << g << std::endl;
                puts("Impossible");
                exit(0);
            }

            c /= g;

            x = x * c ;
            y = y * c ;
        }
    };
}

#endif
