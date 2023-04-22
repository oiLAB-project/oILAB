#ifndef gbLAB_Farey_h_
#define gbLAB_Farey_h_

#include <iostream>
#include <numeric>
#include <vector>

namespace gbLAB {
    // Define a struct to store the current and previous mediant fractions
    struct Fraction {
        int n, d, N, D;

        Fraction(int n_, int d_, int N_, int D_) {
            n = n_;
            d = d_;
            N = N_;
            D = D_;
        }
    };
    std::vector<std::pair<int, int>> farey(int limit, const bool firstQuadrant=true);

}

#endif
