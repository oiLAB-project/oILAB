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
            for(int i=0; i<dim-1; ++i)
            {
                if (this->operator()(i) < rhs(i)) return true;
                if (rhs(i) < this->operator()(i)) return false;
            }
            if (this->operator()(dim-1) < rhs(dim-1)) return true;
            return false;
            /*
            if (this->operator()(0) < rhs(0)) return true;
            if (rhs(0) < this->operator()(0)) return false;
            if (this->operator()(1) < rhs(1)) return true;
            if (rhs(1) < this->operator()(1)) return false;
            if (this->operator()(2) < rhs(2)) return true;
            if (rhs(2) < this->operator()(2)) return false;
            if (this->operator()(3) < rhs(3)) return true;
            return false;
             */
        }

        virtual ~OrderedTuplet() {}
    };

    class XTuplet : public Eigen::Matrix<int,Eigen::Dynamic,1>
    {
    public:
        XTuplet(const int& sz) : Eigen::Matrix<int,Eigen::Dynamic,1>(sz)
        {}

        static void generate_tuples(int n, int k, std::vector<XTuplet>& results, XTuplet& current_tuple, int index) {
            if (index == k) {
                results.push_back(current_tuple);  // Add current tuple to results
            } else {
                for (int i = 0; i < n; ++i) {
                    current_tuple(index) = i;  // Assign current element
                    generate_tuples(n, k, results, current_tuple, index + 1);  // Recursively generate next element
                }
            }
        }

        /*!
         * Constructs the collection of all tuplets of size \p k with each entry in the
         * range \f$\{0,\dots,n-1\}\f$.
         * @param n
         * @param k
         * @param results
         * @param current_tuple
         * @param index
         */
        static std::vector<XTuplet> generate_tuples(int n, int k) {
            std::vector<XTuplet> results;
            XTuplet current_tuple(k);  // Initialize a vector of size k with 0s
            current_tuple.setZero();


            generate_tuples(n, k, results, current_tuple, 0);  // Start generating tuples from index 0

            return results;
        }

        static void generate_tuples(std::vector<int> n, int k, std::vector<XTuplet>& results, XTuplet& current_tuple, int index) {
            if (index == k) {
                results.push_back(current_tuple);  // Add current tuple to results
            } else {
                for (int i = 0; i < n[index]; ++i) {
                    current_tuple(index) = i;  // Assign current element
                    generate_tuples(n, k, results, current_tuple, index + 1);  // Recursively generate next element
                }
            }
        }

        /*!
         * Constructs the collection of all tuplets of size \p k with the i-th entry in the
         * range \f$\{0,\dots,n_i-1\}\f$.
         * @param n
         * @param k
         * @param results
         * @param current_tuple
         * @param index
         */
        static std::vector<XTuplet> generate_tuples(std::vector<int> n, int k) {
            assert(n.size() == k);
            std::vector<XTuplet> results;
            XTuplet current_tuple(k);  // Initialize a vector of size k with 0s
            current_tuple.setZero();


            generate_tuples(n, k, results, current_tuple, 0);  // Start generating tuples from index 0

            return results;
        }

        // Define < operator to use std::map<OrderedTuplet,std::vector<int>>
        bool operator<(const XTuplet& rhs) const {
            for (int i = 0; i < rhs.size()-1; ++i) {
                if (this->operator()(i) < rhs(i)) return true;
                if (rhs(i) < this->operator()(i)) return false;
            }
            if (this->operator()(rhs.size()- 1) < rhs(rhs.size()- 1)) return true;
            return false;
        }


    };

}

#endif //OILAB_ORDEREDTUPLET_H
