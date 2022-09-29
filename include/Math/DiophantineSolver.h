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
//#include <LatticeGCD.h>
#include <IntegerMath.h>
#include<climits>
namespace gbLAB
{
    class DiophantineSolver
    {
        public:
        // Finds such x and y, that a  x + b  y = gcd(a, b)
        //https://github.com/ADJA/algos/blob/master/NumberTheory/DiophantineEquation.cpp
        static long int extended_gcd(long long int a, long long int b, long long int &x, long long int &y)
        {
            if (b == 0)
            {
                x = 1;
                y = 0;
                return a;
            }
            long long int x1, y1;
            long long int g = extended_gcd(b, a % b, x1, y1);
            x = y1;
            y = x1 - (a / b) * y1;
            return g;
        }

        // Solves equation a  x + b  y = c, writes answer to x and y
        static void solveDiophantine2vars(long long int a, long long int b, long long int c, long long int &x, long long int &y)
        {

            long long int g = extended_gcd(a, b, x, y);

            if (c % g != 0)
            {
                std::cout << a << "  " << b << "  " << c << "   " << g << std::endl;
                puts("Impossible");
                exit(0);
            }

            c /= g;
            // long int temp;
            // if (a<b)
            // {
            //     temp=0;
            // }
            // else
            // {
            //     temp=LLONG_MAX;
            // }
            // x = x * c + temp*b/g;
            // y = y * c - temp*a/g;

            x = x * c ;
            y = y * c ;
        }
    };
}

#endif
