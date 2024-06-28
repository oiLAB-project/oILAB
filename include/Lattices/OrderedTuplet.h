//
// Created by Nikhil Chandra Admal on 2/23/24.
//

#ifndef OILAB_ORDEREDTUPLET_H
#define OILAB_ORDEREDTUPLET_H
#include<Eigen/Eigen>

namespace gbLAB {
    template <int dim>
    class OrderedTuplet : public Eigen::Matrix<typename LatticeCore<dim>::IntScalarType,dim,1>
    {
    public:

        OrderedTuplet() = default;

        // Define < operator to use std::map<OrderedTuplet,std::vector<int>>
        bool operator<(const OrderedTuplet &rhs) const {
            if (this->operator()(0) < rhs(0)) return true;
            if (rhs(0) < this->operator()(0)) return false;
            if (this->operator()(1) < rhs(1)) return true;
            if (rhs(1) < this->operator()(1)) return false;
            if (this->operator()(2) < rhs(2)) return true;
            return false;
        }

        static void generate_tuples(int n, vector<vector<int>>& results, vector<int>& current_tuple, int index) {
            if (index == dim) {
                results.push_back(current_tuple);  // Add current tuple to results
            } else {
                for (int i = 1; i <= n; ++i) {
                    current_tuple[index] = i;  // Assign current element
                    generate_tuples(n, results, current_tuple, index + 1);  // Recursively generate next element
                }
            }
        }

        static vector<vector<int>> generate_tuples(int n) {
            vector<vector<int>> results;
            vector<int> current_tuple(dim, 0);  // Initialize a vector of size k with 0s

            generate_tuples(n, results, current_tuple, 0);  // Start generating tuples from index 0

            return results;
        }

        virtual ~OrderedTuplet() {}
    };

}

#endif //OILAB_ORDEREDTUPLET_H
