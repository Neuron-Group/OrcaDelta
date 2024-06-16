
// #include "bits/stdc++.h"
// #include "QUADTREE.hpp"
#include "iostream"
#include "QUADTREE.hpp"
#include "cmath"
#include "RANDOM.hpp"
#include "SOLVER.hpp"

using namespace std;

template<class T>
class myClass {
    public:
    T a;
    T* b;
    myClass(){
        this->b = nullptr;
        return;
    }
    void func(){
        // cout << ISPTR(this->a) << endl;
        return;
    };

    myClass(const myClass &b) {
        this->a = b.a;
    }
    ~myClass(){
        DEL(this->a);
    }
};

int main(){
    size_t T = 50;
    SPARSE::MAT<float, 3> A(T, T);
    A.eye_();
    A.mul_(-2);

    for (size_t i = 0; i < T-1; i++) {
        A.insert(1, i+1, i);
        A.insert(1, i, i+1);
    }

    SPARSE::MAT<float, 3> b(T, 1);
    b.insert(-1, 0, 0);
    b.insert(-1, T-1, 0);
    SPARSE::MAT<float, 3> x(T, 1);

    vector<float> err;

    Gongetidu(A, x, b, err);

    x.test();
}