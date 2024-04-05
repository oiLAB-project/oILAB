//
// Created by Nikhil Chandra Admal on 2/7/24.
//

#ifndef OILAB_RANDOM_H
#define OILAB_RANDOM_H
#include <iostream>
#include <random>

namespace gbLAB {
    template<typename T>
    T random(const T& a, const T& b) {
        std::random_device rd;  // a seed source for the random number engine
        std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd()
        if(typeid(T) == typeid(int)) {
            std::uniform_int_distribution<> distrib(a, b);
            return distrib(gen);
        }
        else if(typeid(T) == typeid(double)) {
            std::uniform_real_distribution<> distrib(a, b);
            return distrib(gen);
        }
        else{
            throw(std::runtime_error("Unknown type\n"));
        }
    }
}
#endif //OILAB_RANDOM_H
