
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

    RandomGenerator T1_rnd(1, 100);
    RandomGenerator T2_rnd(1, 100);

    RandomGenerator ele_rnd(-100, 100);

    int T1 = T1_rnd();
    int T2 = T2_rnd();

    SPARSE::MAT<float, 2> mata(T1, T2, true);
    SPARSE::MAT<float, 2> matb(T2, T1, true);
    SPARSE::MAT<float, 2> matc(T1, T1, true);

    Eigen::MatrixXf mat_a(T1, T2);
    Eigen::MatrixXf mat_b(T2, T1);
    Eigen::MatrixXf mat_c(T1, T1);
    

    for (int i = 0; i < T1; i++) {
        for (int j = 0; j < T2; j ++) {
            int ele_a = ele_rnd();
            mata.insert(ele_a, i, j);
            mat_a(i, j) = ele_a;

            int ele_b = ele_rnd();
            matb.insert(ele_b, j, i);
            mat_b(j, i) = ele_b;
        }
    }

    // 一旦非0层节点中出现节点不满四个孩子节点（孩子节点中存在 nullptr）
    // 乘法就会出错

    // matc = mata;
    matc = mata*matb;
    cout << "1" << endl;
    mat_c = mat_a * mat_b;

    // cout << "mata = " << endl;
    // mata.test();

    // cout << "matb = " << endl;
    // matb.test();

    // cout << "ans = " << endl;
    // matc.test();

    // cout << endl;

    // cout << "ground truth = " << endl;
    // cout << mat_c << endl;




    for (int i = 0; i < T1; i++) {
        for (int j = 0; j < T2; j++) {
            if (mat_c(i, j) != matc.getVal(i, j)) {
                cout << "wrong" << endl;
                return 0;
            }
        }
    }

    cout << "right" << endl;
}