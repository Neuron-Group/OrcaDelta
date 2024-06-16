#include <iostream>
#include <random>

class RandomGenerator {
public:
    RandomGenerator(int l=0, int r=2) : l_(l), r_(r) {
        // 使用系统熵池进行初始随机数生成
        std::random_device rd;
        rng_.seed(rd());
    }

    // 重载括号运算符
    int operator()() {
        // 使用均匀分布生成区间[l, r)的随机数
        std::uniform_int_distribution<int> dist(l_, r_ - 1);
        int rand_num = dist(rng_);
        
        // 更新随机数生成器种子
        rng_.seed(rand_num*234748);
        
        return rand_num;
    }

private:
    int l_;
    int r_;
    std::mt19937 rng_; // 使用Mersenne Twister 19937生成器
};