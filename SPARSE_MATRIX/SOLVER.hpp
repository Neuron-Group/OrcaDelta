
#ifndef __SOLVER_HPP__
#define __SOLVER_HPP__

#pragma GCC optimize (3, "Ofast", "inline")

#include "QUADTREE.hpp"

template<class T, size_t unit_size>
bool Gongetidu(
    SPARSE::MAT<T, unit_size> A,
    SPARSE::MAT<T, unit_size> &x,
    SPARSE::MAT<T, unit_size> b,
    std::vector<T> & error,
    T e0 = 1e-3,
    size_t maxIt = 100
) {
    size_t n = static_cast<size_t>(b.size().getVal(0, 0));
    x = SPARSE::zeros<T, unit_size>(n, 1);
    SPARSE::MAT<T, unit_size> r = b - A*x;
    SPARSE::MAT<T, unit_size> d = r;
    SPARSE::MAT<T, unit_size> alpha, bt;
    SPARSE::MAT<T, unit_size> r1 = r;
    size_t i = 0;

    for (size_t it = 0; it < maxIt; it++) {
        alpha = ((r.transpose())*d)/((d.transpose())*A*d);
        x = x+alpha.getVal(0, 0)*d;
        r1 = b-A*x;
        error.push_back(r1.norm());
        bt = -((r1.transpose())*A*d)/((d.transpose())*A*d);
        d = r1+bt.getVal(0, 0)*d;
        r = r1;
        i = i+1;
        if (r1.norm() <= e0) {
            return true;
        }

        std::cout << it << std::endl;
    }
    return false;
}

#endif //__SOLVER_HPP__
