#include "Farey.h"

namespace gbLAB
{
    // Function to generate the Farey sequence up to a given limit
    std::vector<std::pair<int, int>> farey(int limit, const bool firstQuadrant)
    {
        std::vector<Fraction> pend;
        std::vector<std::pair<int, int>> output;
        int n = 0;
        int d = 1;
        int N = 1;
        int D = 1;
        while (true) {
            int mediant_d = d + D;
            if (mediant_d <= limit) {
                int mediant_n = n + N;
                pend.emplace_back(Fraction(mediant_n, mediant_d, N, D));
                N = mediant_n;
                D = mediant_d;
            } else {
                output.emplace_back(n, d);
                if (firstQuadrant == false) {
                    if (n != 0) output.emplace_back(-n, d);
                    if (d != 0) output.emplace_back(n, -d);
                    if (n != 0 && d != 0) output.emplace_back(-n, -d);
                    output.emplace_back(d, n);
                    if (n != 0) output.emplace_back(d, -n);
                    if (d != 0) output.emplace_back(-d, n);
                    if (n != 0 && d != 0) output.emplace_back(-d, -n);
                }
                if (pend.empty()) {
                    break;
                }
                Fraction frac = pend.back();
                pend.pop_back();
                n = frac.n;
                d = frac.d;
                N = frac.N;
                D = frac.D;
            }
        }
        output.emplace_back(1, 1);
        return output;
    }
}

