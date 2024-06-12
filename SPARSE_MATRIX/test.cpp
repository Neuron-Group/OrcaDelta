
// #include "bits/stdc++.h"
// #include "QUADTREE.hpp"
#include "iostream"
#include "QUADTREE.hpp"
#include "cmath"

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
    int T = 51;

    SPARSE::MAT<float, 3> mata(T, T, true);

    SPARSE::MAT<float, 3> matb(T, T, true);
    matb = mata.transpose();


    matb = mata.setZero();
    mata.setZero_();

    for (int i = 0; i < T; i++) {
        mata.insert(i, i, i);
    }
    for (int j = 0; j < T; j++) {
        matb.insert(j, T-j-1, j);
    }

    cout << "mata = " << endl;
    mata.test();

    cout << "matb = " << endl;
    matb.test();

    matb = matb - mata;
    mata = matb + mata;

    cout << "matb - mata + mata = matb = " << endl;
    mata.test();

    cout << "matb - mata = " << endl;
    matb.test();

    mata.neg_();

    cout << "-matb = " << endl;
    mata.test();

    matb = -mata;

    cout << "matb = " << endl;
    matb.test();
    // cout << endl;
    /*
    SPARSE::MAT<float, 3> vec(10, 1);
    for (int i = 0; i < 10; i++) {
        vec.insert(i, i, 0);
    }
    for (int i = 0; i < 10; i++) {
        // cout << vec.getVal(i, 0) << endl;
    }
    // cout << endl;

    for (int i = 0; i < 10; i++) {
        // cout << vec.transpose().getVal(0, i) << " ";
    }
    // cout << endl;
    */

    cout << mata.size().getVal(0, 0) << endl;
    cout << mata.size().getVal(1, 0) << endl;
}