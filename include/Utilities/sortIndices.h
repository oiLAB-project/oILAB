//
// Created by Nikhil Chandra Admal on 9/30/22.
// https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes

#ifndef OILAB_SORTINDICES_H
#define OILAB_SORTINDICES_H

#include <iostream>
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort

using namespace std;

namespace gbLAB {
    template<typename ArrayType>
    vector <size_t> sortIndices(const ArrayType &v) {

        // initialize original index locations
        vector<size_t> idx(v.size());
        iota(idx.begin(), idx.end(), 0);

        // sort indexes based on comparing values in v
        // using std::stable_sort instead of std::sort
        // to avoid unnecessary index re-orderings
        // when v contains elements of equal values
        stable_sort(idx.begin(), idx.end(),
                    [&v](size_t i1, size_t i2) { return abs(v(i1)) > abs(v(i2)); });

        return idx;
    }
}
#endif //OILAB_SORTINDICES_H
