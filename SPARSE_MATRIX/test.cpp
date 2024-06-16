
// #include "bits/stdc++.h"
// #include "QUADTREE.hpp"
#include "iostream"
#include "QUADTREE.hpp"
#include "cmath"
#include "RANDOM.hpp"

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

    RandomGenerator T1_rnd(1, 10);
    RandomGenerator T2_rnd(1, 10);

    RandomGenerator ele_rnd(-50, 50);

    int T1 = T1_rnd();
    int T2 = T2_rnd();

    SPARSE::MAT<float, 3> mata(T1, T2, true);
    SPARSE::MAT<float, 3> matb(T1, T2, true);
    SPARSE::MAT<float, 3> matc(T1, T1, true);

    Eigen::MatrixXf mat_a(T1, T2);
    Eigen::MatrixXf mat_b(T2, T1);
    Eigen::MatrixXf mat_c(T1, T1);
    

    for (int i = 0; i < T1; i++) {
        for (int j = 0; j < T2; j ++) {
            int ele_a = ele_rnd();
            mata.insert(static_cast<float>(i+1), i, j);
            // mat_a(i, j) = static_cast<float>(ele_a);

            int ele_b = ele_rnd();
            matb.insert(static_cast<float>(j+1), i, j);
            // mat_b(j, i) = static_cast<float>(ele_b);
        }
    }
    mata.test();
    matb.test();
    SPARSE::divWithEle(mata, matb).test();

    SPARSE::eye<float, 3>(7).test();
}