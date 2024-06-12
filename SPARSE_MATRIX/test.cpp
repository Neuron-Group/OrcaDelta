
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
        cout << ISPTR(this->a) << endl;
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

    SPARSE::MAT<float, 3> mata(T, T);

    for (int i = 0; i < T; i++) {
        for (int j = 0; j < T; j ++){
            mata.insert(i*T + j, i, j);
        }
    }

    for (int i = 0; i < T; i++) {
        for (int j = 0; j < T; j++){
            cout << mata.getVal(i, j) << " ";
        }
        cout << endl;
    }
    cout << endl;
    SPARSE::MAT<float, 3> matb(T, T);
    matb = mata.transpose();

    for (int i = 0; i < T; i++) {
        for (int j = 0; j < T; j++){
            cout << matb.getVal(i, j) << " ";
        }
        cout << endl;
    }

    cout << endl;
    for (int i = 0; i < T; i++) {
        for (int j = 0; j < T; j++){
            cout << (mata + matb).getVal(i, j) << " ";
        }
        cout << endl;
    }

    matb = mata.setZero();
    mata.setZero_();

    for (int i = 0; i < T; i++) {
        for (int j = 0; j < T; j ++){
            if (i == j) mata.insert(i, i, j);
        }
    }
    for (int i = 0; i < T; i++) {
        for (int j = 0; j < T; j ++){
            if (i == T-j-1) matb.insert(i, i, j);
        }
    }

    for (int i = 0; i < T; i++) {
        for (int j = 0; j < T; j++){
            cout << mata.getVal(i, j) << " ";
        }
        cout << endl;
    }
    cout << endl;

    for (int i = 0; i < T; i++) {
        for (int j = 0; j < T; j++){
            cout << matb.getVal(i, j) << " ";
        }
        cout << endl;
    }

    cout << endl;
    for (int i = 0; i < T; i++) {
        for (int j = 0; j < T; j++){
            cout << (mata - matb).getVal(i, j) << " ";
        }
        cout << endl;
    }
    /*
    SPARSE::MAT<float, 3> vec(10, 1);
    for (int i = 0; i < 10; i++) {
        vec.insert(i, i, 0);
    }
    for (int i = 0; i < 10; i++) {
        cout << vec.getVal(i, 0) << endl;
    }
    cout << endl;

    for (int i = 0; i < 10; i++) {
        cout << vec.transpose().getVal(0, i) << " ";
    }
    cout << endl;
    */
}