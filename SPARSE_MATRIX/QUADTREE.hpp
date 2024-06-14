
#ifndef __QUADTREE_HPP__
#define __QUADTREE_HPP__

#include "iostream"
#include "cstdio"
#include "algorithm"
#include "cmath"

#include "queue"
#include "vector"
#include "stack"

#include "eigen3/Eigen/Dense"

template<typename T>
char IsPtr(T* v)
{
    return 'd';
}
int IsPtr(...)				// 三个点表示可接受任何函数
{
    return 0;
}
#define ISPTR(p) (sizeof(IsPtr(p)) == sizeof(char))

template<typename T>
int DEL(T* v)
{
    delete v;
    return 1;
}
int DEL(...)				// 三个点表示可接受任何函数
{
    return 0;
}

namespace SPARSE{

    /*
        ### 示例：自定义异常类

        下面是一个示例，演示如何创建和使用自定义的异常类：

        ```cpp
        #include <iostream>
        #include <exception>
        #include <stdexcept>
        #include <string>

        // 自定义异常类继承自 std::runtime_error
        class IndexOutOfRangeException : public std::runtime_error {
        public:
            explicit IndexOutOfRangeException(const std::string& message)
                : std::runtime_error(message) {}
        };

        void functionThatThrows() {
            throw IndexOutOfRangeException("Index out of range.");
        }

        int main() {
            try {
                functionThatThrows();
            } catch (const IndexOutOfRangeException& e) {
                std::cerr << "Caught an exception: " << e.what() << std::endl;
            } catch (const std::exception& e) {
                // 捕捉其他继承自 std::exception 的异常
                std::cerr << "Caught an std::exception: " << e.what() << std::endl;
            } catch (...) {
                // 捕捉所有其他类型的异常
                std::cerr << "Caught an unknown exception." << std::endl;
            }

            return 0;
        }
        ```

        ### 解释：
        1. **自定义异常类**：`IndexOutOfRangeException` 继承自 `std::runtime_error`，这是一个标准库提供的异常类，专门用于运行时错误。
        
        2. **抛出异常**：在 `functionThatThrows` 函数中抛出 `IndexOutOfRangeException`，并传递一个错误消息。

        3. **捕捉异常**：在 `main` 函数中使用 `try-catch` 语句捕捉异常，并使用 `e.what()` 来获得错误消息的内容。

        当自定义异常类并使用继承自 `std::exception` 的异常类时，`what()` 方法会返回你在抛出异常时传递的错误消息，从而可以更加清晰地描述出现的问题。这样一来，你在捕捉和处理异常时就能获得更详细的信息，而不仅仅是类型信息。

        通过这种方式，当你执行程序并抛出 `"Index out of range."` 异常时，你会看到更有意义的错误信息，如：

        ```
        Caught an exception: Index out of range.
        ```
    */

    // 下标溢出异常
    // 自定义异常类继承自 std::runtime_error
    class IndexOutOfRangeException : public std::runtime_error {
        public:
        explicit IndexOutOfRangeException(const std::string& message)
            : std::runtime_error(message) {}
    };

    void functionThatThrows() {
        throw IndexOutOfRangeException("Index out of range.");
    }

    // 树节点类
    template <class T>
    class TREENODE {

        size_t deepth;  // 从下往上的深度
        T * dataNode;
        TREENODE * CHILD[4];

        public:
        
        template<typename U, size_t unit_size>
        friend class MAT;

        TREENODE(size_t deepth = 0);
        TREENODE<T> & operator = (const TREENODE<T> &b);
        TREENODE(const TREENODE &b);
        ~TREENODE() {
            if (this->dataNode == nullptr) return;
            delete[] this->dataNode;
        }

        void insert(TREENODE<T> * & root, size_t deepth, T DATA, size_t i, size_t j){
            if(root == nullptr){
                root = new TREENODE<T>(deepth);
            }
            TREENODE * ptr = root;
            while (ptr->deepth > 0) {
                short k = 0;
                if (i & (1<<(ptr->deepth-1))){
                    k = k | 1;
                }
                if (j & (1<<(ptr->deepth-1))){
                    k = k | 2;
                }
                if (ptr->CHILD[k] == nullptr){
                    ptr->CHILD[k] = new TREENODE(ptr->deepth - 1);
                }
                ptr = ptr->CHILD[k];
            }
            if (ptr->dataNode == nullptr){
                ptr->dataNode = new T(DATA);
            } else {
                *(ptr->dataNode) = DATA;
            }
        }

        T getVal(){
            return *(this->dataNode);
        }

        TREENODE * get(size_t i, size_t j){
            TREENODE * ptr = this;
            if(ptr == nullptr){
                return ptr;
            }
            while (ptr->deepth > 0) {
                short k = 0;
                if (i & (1<<(ptr->deepth-1))){
                    k = k | 1;
                }
                if (j & (1<<(ptr->deepth-1))){
                    k = k | 2;
                }
                if (ptr->CHILD[k] == nullptr){
                    return nullptr;
                }
                ptr = ptr->CHILD[k];
            }
            return ptr;
        }

        void test_print(TREENODE<T>* root, short k_ = 0, size_t level = 0) {
            if (root == nullptr) {
                return;
            }
            for (long long i = 0; i < static_cast<long long>(level-1); i++) {
                std::cout << "│  ";
            }
            if (level > 0) std::cout << "├──";
            std::cout << "deepth : " << root->deepth << ", position : " << k_ << std::endl;
            for (short i = 0; i < 4; i++) {
                test_print(root->CHILD[i], i, level+1);
            } 
        }

        void add_(const TREENODE<T> &b);
        TREENODE<T> operator + (const TREENODE<T> &b) const;

        void sub_(const TREENODE<T> &b);
        TREENODE<T> operator - (const TREENODE<T> &b) const;

        void neg_();
        TREENODE<T> operator - () const;
    };

    template <class T>
    TREENODE<T>::TREENODE(size_t deepth){
        this->deepth = deepth;
        this->dataNode = nullptr;
        this->CHILD[0] = nullptr;
        this->CHILD[1] = nullptr;
        this->CHILD[2] = nullptr;
        this->CHILD[3] = nullptr;
    }

    template<class T>
    TREENODE<T> & TREENODE<T>::operator = (const TREENODE<T> &b) {

        if (this == &b) return *this;

        if (b.dataNode == nullptr) {
            delete this->dataNode;
            this->dataNode = nullptr;
        } else {
            if (this->dataNode == nullptr) {
                this->dataNode = new T(*(b.dataNode));
            } else {
                *(this->dataNode) = *(b.dataNode);
            }
        }

        this->deepth = b.deepth;

        return *this;
    }

    template<class T>
    TREENODE<T>::TREENODE(const TREENODE &b) {
        this->deepth = b.deepth;
        this->dataNode = nullptr;
        this->CHILD[0] = nullptr;
        this->CHILD[1] = nullptr;
        this->CHILD[2] = nullptr;
        this->CHILD[3] = nullptr;

        if (b.dataNode == nullptr) {
            delete this->dataNode;
            this->dataNode = nullptr;
        } else {
            this->dataNode = new T(*(b.dataNode));
        }
    }

    template<class T>
    void TREENODE<T>::add_(const TREENODE<T> & b) {
        if (b.dataNode == nullptr) return;
        if (this->dataNode == nullptr) {
            this->dataNode = new T(*(b.dataNode));
            return;
        }
        *(this->dataNode) = *(this->dataNode) + *(b.dataNode);
    }

    template<class T>
    TREENODE<T> TREENODE<T>::operator + (const TREENODE<T> & b) const {
        TREENODE<T> ans = *this;
        ans.add_(b);
        return ans;
    }

    template<class T>
    void TREENODE<T>::sub_(const TREENODE<T> & b) {
        if (b.dataNode == nullptr) return;
        if (this->dataNode == nullptr) {
            this->dataNode = new T(-*(b.dataNode));
            return;
        }
        *(this->dataNode) = *(this->dataNode) - *(b.dataNode);
    }

    template<class T>
    TREENODE<T> TREENODE<T>::operator - (const TREENODE<T> & b) const {
        TREENODE<T> ans = *this;
        ans.add_(b);
        return ans;
    }

    template<class T>
    void TREENODE<T>::neg_() {
        if (this->dataNode == nullptr) return;
        *(this->dataNode) = -*(this->dataNode);
    }

    template<class T>
    TREENODE<T> TREENODE<T>::operator - () const {
        TREENODE<T> ans = *this;
        ans.neg_();
        return ans;
    }

    // 矩阵类
    template<class T, size_t unit_size>
    class MAT {
        size_t rows, cols;

        bool autoformat; // 是否自动删除 0 元素
        long double prec; // 判 0 精度

        TREENODE<Eigen::Matrix<T, unit_size, unit_size>> * ROOT;

        // 判 0 函数
        // 矩阵判0
        bool iszero_mat(const Eigen::Matrix<T, unit_size, unit_size> & mat) const {
            for (size_t i = 0; i < unit_size; i++) {
                for (size_t j = 0; j < unit_size; j++) {
                    if (static_cast<long double>(abs(mat(i, j))) > this->prec) {
                        return false;
                    }
                }
            }
            return true;
        }

        // 标量判 0
        bool iszero_sca(const T x) const {
            if (static_cast<long double>(abs(x)) > this->prec) return false;
            return true;
        }

        size_t quadtree_size() {
            return (std::max(this->rows, this->cols) + unit_size - 1)/unit_size;
        }

        size_t quadtree_high() {
            return std::ceil(log2(quadtree_size()))+1;
        }

        size_t IDX(size_t idx) {
            return idx/unit_size;
        }
        
        size_t IDX_(size_t idx) {
            return idx%unit_size;
        }

        bool check_idx(size_t i, size_t j){
            return (i < this->rows && j < this->cols) ? true : false;
        }

        void check_idx_with_error(size_t i, size_t j) {
            if (check_idx(i, j)) return;
            throw IndexOutOfRangeException("Index out of range.");
        }

        void initForData(TREENODE<Eigen::Matrix<T, unit_size, unit_size>> * node) {
            if (node->dataNode == nullptr) {
                node->dataNode = new Eigen::Matrix<T, unit_size, unit_size>;
                node->dataNode->setZero();
            }
        }

        // 移除子树
        void del(TREENODE<Eigen::Matrix<T, unit_size, unit_size>> * & root) {
            if (!(root == nullptr)) {
                std::queue<TREENODE<Eigen::Matrix<T, unit_size, unit_size>> *> q;
                q.push(root);
                while (!q.empty()) {
                    TREENODE<Eigen::Matrix<T, unit_size, unit_size>> * f = q.front();
                    if (f->CHILD[0] != nullptr) q.push(f->CHILD[0]);
                    if (f->CHILD[1] != nullptr) q.push(f->CHILD[1]);
                    if (f->CHILD[2] != nullptr) q.push(f->CHILD[2]);
                    if (f->CHILD[3] != nullptr) q.push(f->CHILD[3]);
                    // delete f->dataNode;
                    delete f;
                    q.pop();
                }
                root= nullptr;
            }
        }

        // 移除单个节点
        void del_single_node(TREENODE<Eigen::Matrix<T, unit_size, unit_size>> * & node) {
            delete node;
            node = nullptr;
        }

        // 删除 0 元素及其分支
        short removezero(TREENODE<Eigen::Matrix<T, unit_size, unit_size>> * & root) {
            if (root == nullptr) {
                return 0;
            }

            if (root->deepth == 0) {
                if (root->dataNode == nullptr) {
                    del_single_node(root);
                    return 0;
                } else {
                    if (iszero_mat(*(root->dataNode))) {
                        del_single_node(root);
                        return 0;
                    } else {
                        return 1;
                    }
                }
            }

            short flag = 0;

            for (short i = 0; i < 4; i++) {
                flag += removezero(root->CHILD[i]);
            }

            if (flag == 0) {
                del_single_node(root);
                return 0;
            } else {
                return 1;
            }
        }

        public:
        MAT(size_t rows = 1, size_t cols = 1, bool autoformat = true, long double prec = 1e-8);
        MAT<T, unit_size> & operator = (const MAT<T, unit_size> &b);
        MAT(const MAT &b);
        ~MAT();

        bool isBlank () const {
            return (this->ROOT == nullptr) ? true : false;
        }

        /*
        bool insert(T data, size_t i, size_t j){
            check_idx_with_error(i, j);
            // if (this->iszero_sca(data)) return true;
            if (check_idx(i, j) != true) return false;
            size_t I = IDX(i);
            size_t J = IDX(j);
            if(isBlank()){
                this->ROOT = new TREENODE<Eigen::Matrix<T, unit_size, unit_size>>(this->quadtree_high()-1);
                if (this->ROOT->deepth == 0) {
                    this->initForData(ROOT);
                }
            }
            TREENODE<Eigen::Matrix<T, unit_size, unit_size>> * ptr = ROOT;
            while (ptr->deepth > 0) {
                short k = 0;
                if (I & (1<<(ptr->deepth-1))){
                    k = k | 1;
                }
                if (J & (1<<(ptr->deepth-1))){
                    k = k | 2;
                }
                if (ptr->CHILD[k] == nullptr){
                    ptr->CHILD[k] = new TREENODE<Eigen::Matrix<T, unit_size, unit_size>>(ptr->deepth - 1);
                }
                ptr = ptr->CHILD[k];
            }
            this->initForData(ptr);
            (*(ptr->dataNode))(IDX_(i), IDX_(j)) = data;
            return true;
        }
        */

       bool insert(T data, size_t i, size_t j){
            check_idx_with_error(i, j);
            // if (this->iszero_sca(data)) return true;
            if (check_idx(i, j) != true) return false;
            size_t I = IDX(i);
            size_t J = IDX(j);
            if(isBlank()){
                this->ROOT = new TREENODE<Eigen::Matrix<T, unit_size, unit_size>>(this->quadtree_high()-1);
                if (this->ROOT->deepth == 0) {
                    this->initForData(ROOT);
                }
            }
            std::stack<TREENODE<Eigen::Matrix<T, unit_size, unit_size>> *> ptr_stk;
            std::stack<short> dir_stk;

            TREENODE<Eigen::Matrix<T, unit_size, unit_size>> * ptr = ROOT;
            // ptr_stk.push(ROOT);
            while (ptr->deepth > 0) {
                short k = 0;
                if (I & (1<<(ptr->deepth-1))){
                    k = k | 1;
                }
                if (J & (1<<(ptr->deepth-1))){
                    k = k | 2;
                }
                if (ptr->CHILD[k] == nullptr){
                    ptr->CHILD[k] = new TREENODE<Eigen::Matrix<T, unit_size, unit_size>>(ptr->deepth - 1);
                }
                ptr_stk.push(ptr->CHILD[k]);
                dir_stk.push(k);
                ptr = ptr->CHILD[k];
            }
            this->initForData(ptr);
            (*(ptr->dataNode))(IDX_(i), IDX_(j)) = data;

            // 清除 0 元素
            if (iszero_mat(*(ptr->dataNode)) && this->autoformat == true) {

                // ptr = ptr_stk.top();
                if (ptr == ROOT) {
                    delete ptr;
                    ROOT = nullptr;
                    return true;
                }
                delete ptr;
                ptr_stk.pop();

                short k;

                while(!ptr_stk.empty()) {
                    
                    ptr = ptr_stk.top();
                    k = dir_stk.top();
                    ptr->CHILD[k] = nullptr;
                    
                    short child_cnt = 0;
                    if(ptr->CHILD[0] != nullptr) child_cnt++;
                    if(ptr->CHILD[1] != nullptr) child_cnt++;
                    if(ptr->CHILD[2] != nullptr) child_cnt++;
                    if(ptr->CHILD[3] != nullptr) child_cnt++;
                    
                    if (child_cnt != 0) {
                        return true;
                    }

                    delete ptr;

                    ptr_stk.pop();
                    dir_stk.pop();
                }

                k = dir_stk.top();
                ROOT->CHILD[k] = nullptr;

                short child_cnt = 0;
                if(ROOT->CHILD[0] != nullptr) child_cnt++;
                if(ROOT->CHILD[1] != nullptr) child_cnt++;
                if(ROOT->CHILD[2] != nullptr) child_cnt++;
                if(ROOT->CHILD[3] != nullptr) child_cnt++;

                if (child_cnt == 0) {
                    del_single_node(this->ROOT);
                }
            }
            return true;
        }

        T getVal(size_t i, size_t j) {
            check_idx_with_error(i, j);
            // 索引到块
            TREENODE<Eigen::Matrix<T, unit_size, unit_size>> * ptr = this->ROOT->get(IDX(i), IDX(j));
            
            // 在块中索引出元素
            if (ptr == nullptr) return static_cast<T>(0);
            if (ptr->dataNode == nullptr) return static_cast<T>(0);
            return (*(ptr->dataNode))(IDX_(i), IDX_(j));
        }

        // 格式化
        void format() {
            this->removezero(this->ROOT);
        }

        void test() {
            std::cout << std::endl;
            ROOT->test_print(ROOT);
            // cout << endl;
            for (size_t i = 0; i < this->rows; i++) {
                for (size_t j = 0; j < this->cols; j++) {
                    std::cout << getVal(i, j) << " ";
                }
                std::cout << std::endl;
            }
        }

        // 维度输出
        MAT<long long int, unit_size> size() const {
            MAT<long long int, unit_size> ans(2, 1);
            ans.insert(this->rows, 0, 0);
            ans.insert(this->cols, 1, 0);
            return ans;
        }

        // 打印精度
        long double getPrec() const {
            return this->prec;
        }

        // 设置精度
        void setPrec(const long double prec) {this->prec = prec;}

        // 是否自动稀疏化
        bool isAutoFormat() const {
            return this->autoformat;
        }

        // 设置自动稀疏化
        void setAutoFormat(const bool autoFormat) {
            this->autoformat = autoFormat;
        }

        // 置零
        void setZero_();
        MAT<T, unit_size> setZero() const;

        // 转置
        void transpose_();
        MAT<T, unit_size> transpose() const;

        // 矩阵按位和
        void add_(const MAT<T, unit_size> &b);
        MAT<T, unit_size> operator + (const MAT<T, unit_size> &b) const;

        // 矩阵按位差
        void sub_(const MAT<T, unit_size> &b);
        MAT<T, unit_size> operator - (const MAT<T, unit_size> &b) const;

        // 重载负号
        void neg_();
        MAT<T, unit_size> operator - () const;
    };

    template<class T, size_t unit_size>
    MAT<T, unit_size>::MAT(size_t rows, size_t cols, bool autoformat, long double prec) : 
    rows(rows), 
    cols(cols),
    prec(prec),
    autoformat(autoformat)
    {
        this->ROOT = nullptr;
    }

    template<class T, size_t unit_size>
    MAT<T, unit_size>::~MAT(){
        if (!(this->isBlank())) {
            std::queue<TREENODE<Eigen::Matrix<T, unit_size, unit_size>> *> q;
            q.push(this->ROOT);
            while (!q.empty()) {
                TREENODE<Eigen::Matrix<T, unit_size, unit_size>> * f = q.front();
                if (f->CHILD[0] != nullptr) q.push(f->CHILD[0]);
                if (f->CHILD[1] != nullptr) q.push(f->CHILD[1]);
                if (f->CHILD[2] != nullptr) q.push(f->CHILD[2]);
                if (f->CHILD[3] != nullptr) q.push(f->CHILD[3]);
                // delete f->dataNode;
                delete f;
                q.pop();
            }
            this->ROOT = nullptr;
        }
        
    }

    template<class T, size_t unit_size>
    MAT<T, unit_size> & MAT<T, unit_size>::operator = (const MAT<T, unit_size> &b) {

        this->prec = b.prec;
        this->autoformat = b.autoformat;

        if (this == &b) return *this;

        this->cols = b.cols;
        this->rows = b.rows;
        if (b.isBlank()) {
            del(this->ROOT);
            return *(this);
        }
        if (this->isBlank()) {
            this->ROOT = new TREENODE<Eigen::Matrix<T, unit_size, unit_size>>(b.ROOT->deepth);
            *(this->ROOT) = *(b.ROOT);
        }

        std::queue<TREENODE<Eigen::Matrix<T, unit_size, unit_size>> *> q_ptr_this;
        std::queue<TREENODE<Eigen::Matrix<T, unit_size, unit_size>> *> q_ptr_b;
        
        TREENODE<Eigen::Matrix<T, unit_size, unit_size>> * ptr_this;
        TREENODE<Eigen::Matrix<T, unit_size, unit_size>> * ptr_b;
        
        q_ptr_this.push(this->ROOT);
        q_ptr_b.push(b.ROOT);

        while(!q_ptr_this.empty()) {
            ptr_this = q_ptr_this.front();
            ptr_b = q_ptr_b.front();

            if (ptr_this->deepth == 0) {
                (*(ptr_this)) = (*(ptr_b));
            } else {
                for (short i = 0; i < 4; i++){
                    if (ptr_b->CHILD[i] != nullptr) {
                        if (ptr_this->CHILD[i] == nullptr) {
                            ptr_this->CHILD[i] = new TREENODE<Eigen::Matrix<T, unit_size, unit_size>>(
                                    ptr_this->deepth - 1
                                );
                        }
                        q_ptr_this.push(ptr_this->CHILD[i]);
                        q_ptr_b.push(ptr_b->CHILD[i]);
                    } else {
                        del(ptr_this->CHILD[i]);
                    }
                }
            }

            q_ptr_this.pop();
            q_ptr_b.pop();
        }

        return *this;
    }

    template<class T, size_t unit_size>
    MAT<T, unit_size>::MAT(const MAT &b) {
        this->rows = b.rows;
        this->cols = b.cols;
        this->ROOT = nullptr;
        this->prec = b.prec;
        this->autoformat = b.autoformat;
        if (b.isBlank()) {
            del(this->ROOT);
            return;
        }
        if (this->isBlank()) {
            this->ROOT = new TREENODE<Eigen::Matrix<T, unit_size, unit_size>>(b.ROOT->deepth);
            *(this->ROOT) = *(b.ROOT);
        }

        std::queue<TREENODE<Eigen::Matrix<T, unit_size, unit_size>> *> q_ptr_this;
        std::queue<TREENODE<Eigen::Matrix<T, unit_size, unit_size>> *> q_ptr_b;

        TREENODE<Eigen::Matrix<T, unit_size, unit_size>> * ptr_this;
        TREENODE<Eigen::Matrix<T, unit_size, unit_size>> * ptr_b;
        
        q_ptr_this.push(this->ROOT);
        q_ptr_b.push(b.ROOT);
        while(!q_ptr_this.empty()) {
            ptr_this = q_ptr_this.front();
            ptr_b = q_ptr_b.front();

            if (ptr_this->deepth == 0) {
                (*(ptr_this)) = (*(ptr_b));
            } else {
                for (short i = 0; i < 4; i++){
                    if (ptr_b->CHILD[i] != nullptr) {
                        ptr_this->CHILD[i] = new TREENODE<Eigen::Matrix<T, unit_size, unit_size>>(
                                ptr_this->deepth - 1
                            );
                        q_ptr_this.push(ptr_this->CHILD[i]);
                        q_ptr_b.push(ptr_b->CHILD[i]);
                    }
                }
            }

            q_ptr_this.pop();
            q_ptr_b.pop();
        }
    }

    template<class T, size_t unit_size>
    void MAT<T, unit_size>::setZero_() {
        del(this->ROOT);
        this->ROOT = nullptr;
    }

    template<class T, size_t unit_size>
    MAT<T, unit_size> MAT<T, unit_size>::setZero() const {
        MAT<T, unit_size> ans(this->rows, this->cols);
        return ans;
    }

    template<class T, size_t unit_size>
    void MAT<T, unit_size>::transpose_() {

        size_t tmp_size = this->cols;
        this->cols = this->rows;
        this->rows = tmp_size;

        if (isBlank()) return;

        std::queue<TREENODE<Eigen::Matrix<T, unit_size, unit_size>> *> q;
        q.push(this->ROOT);
        TREENODE<Eigen::Matrix<T, unit_size, unit_size>> * f;
        TREENODE<Eigen::Matrix<T, unit_size, unit_size>> * tmp;
        while (!q.empty()) {
            f = q.front();

            if (f->deepth == 0) {
                if (f->dataNode != nullptr) f->dataNode->transposeInPlace();
            } else {
                tmp = f->CHILD[1];
                f->CHILD[1] = f->CHILD[2];
                f->CHILD[2] = tmp;

                if (f->CHILD[0] != nullptr) q.push(f->CHILD[0]);
                if (f->CHILD[1] != nullptr) q.push(f->CHILD[1]);
                if (f->CHILD[2] != nullptr) q.push(f->CHILD[2]);
                if (f->CHILD[3] != nullptr) q.push(f->CHILD[3]);
            }

            q.pop();
        }
    }

    template<class T, size_t unit_size>
    MAT<T, unit_size> MAT<T, unit_size>::transpose() const {
        MAT<T, unit_size> ans = *this;
        ans.transpose_();
        return ans;
    }

    template<class T, size_t unit_size>
    void MAT<T, unit_size>::add_(const MAT<T, unit_size> & b) {
        // check for size
        if (this->cols != b.cols || this->rows != b.rows) {
            IndexOutOfRangeException("Wrong size!");
        }

        if (b.isBlank()) return;

        if (this->isBlank()) {
            this->ROOT = new TREENODE<Eigen::Matrix<T, unit_size, unit_size>>(*(b.ROOT));
        }

        std::queue<TREENODE<Eigen::Matrix<T, unit_size, unit_size>> *> q_ptr_this;
        std::queue<TREENODE<Eigen::Matrix<T, unit_size, unit_size>> *> q_ptr_b;
        
        TREENODE<Eigen::Matrix<T, unit_size, unit_size>> * ptr_this;
        TREENODE<Eigen::Matrix<T, unit_size, unit_size>> * ptr_b;
        
        q_ptr_this.push(this->ROOT);
        q_ptr_b.push(b.ROOT);
        
        while(!q_ptr_this.empty()) {

            ptr_this = q_ptr_this.front();
            ptr_b = q_ptr_b.front();

            if (ptr_this->deepth == 0) {
                (*(ptr_this)).add_(*(ptr_b));
            } else {
                for (short i = 0; i < 4; i++){
                    if(ptr_b->CHILD[i] != nullptr) {
                        if (ptr_this->CHILD[i] == nullptr) {
                            ptr_this->CHILD[i] = new TREENODE<Eigen::Matrix<T, unit_size, unit_size>>(
                                ptr_this->deepth - 1
                            );
                        }
                        q_ptr_this.push(ptr_this->CHILD[i]);
                        q_ptr_b.push(ptr_b->CHILD[i]);
                    }
                }
            }

            q_ptr_this.pop();
            q_ptr_b.pop();
        }
        if (autoformat) format();
    }

    template<class T, size_t unit_size>
    MAT<T, unit_size> MAT<T, unit_size>::operator+ (const MAT<T, unit_size> & b) const {
        MAT<T, unit_size> ans = *this;
        ans.add_(b);
        return ans;
    }

    template<class T, size_t unit_size>
    void MAT<T, unit_size>::sub_(const MAT<T, unit_size> & b) {

        // check for size
        if (this->cols != b.cols || this->rows != b.rows) {
            IndexOutOfRangeException("Wrong size!");
        }

        if (b.isBlank()) return;

        if (this->isBlank()) {
            this->ROOT = new TREENODE<Eigen::Matrix<T, unit_size, unit_size>>(b.ROOT->deepth);
            *(this->ROOT) = -*(b.ROOT);
        }

        std::queue<TREENODE<Eigen::Matrix<T, unit_size, unit_size>> *> q_ptr_this;
        std::queue<TREENODE<Eigen::Matrix<T, unit_size, unit_size>> *> q_ptr_b;
        
        TREENODE<Eigen::Matrix<T, unit_size, unit_size>> * ptr_this;
        TREENODE<Eigen::Matrix<T, unit_size, unit_size>> * ptr_b;
        
        q_ptr_this.push(this->ROOT);
        q_ptr_b.push(b.ROOT);
        
        while(!q_ptr_this.empty()) {

            ptr_this = q_ptr_this.front();
            ptr_b = q_ptr_b.front();

            if (ptr_this->deepth == 0) {
                (*(ptr_this)).sub_(*(ptr_b));
            } else {
                for (short i = 0; i < 4; i++){
                    if(ptr_b->CHILD[i] != nullptr) {
                        if (ptr_this->CHILD[i] == nullptr) {
                            ptr_this->CHILD[i] = new TREENODE<Eigen::Matrix<T, unit_size, unit_size>>(
                                ptr_this->deepth - 1
                            );
                        }
                        q_ptr_this.push(ptr_this->CHILD[i]);
                        q_ptr_b.push(ptr_b->CHILD[i]);
                    }
                }
            }

            q_ptr_this.pop();
            q_ptr_b.pop();
        }
        if (autoformat) format();
    }

    template<class T, size_t unit_size>
    MAT<T, unit_size> MAT<T, unit_size>::operator- (const MAT<T, unit_size> & b) const {
        MAT<T, unit_size> ans = *this;
        ans.sub_(b);
        return ans;
    }

    template<class T, size_t unit_size>
    void MAT<T, unit_size>::neg_() {
        if (isBlank()) return;
        std::queue<TREENODE<Eigen::Matrix<T, unit_size, unit_size>> *> q;
        q.push(this->ROOT);

        TREENODE<Eigen::Matrix<T, unit_size, unit_size>> * ptr;

        while(!q.empty()) {
            ptr = q.front();

            if (ptr->deepth == 0) {
                if (ptr->dataNode != nullptr) {
                    ptr->neg_();
                }
            } else {
                if (ptr->CHILD[0] != nullptr) q.push(ptr->CHILD[0]);
                if (ptr->CHILD[1] != nullptr) q.push(ptr->CHILD[1]);
                if (ptr->CHILD[2] != nullptr) q.push(ptr->CHILD[2]);
                if (ptr->CHILD[3] != nullptr) q.push(ptr->CHILD[3]);
            }

            q.pop();
        }
    }

    template<class T, size_t unit_size>
    MAT<T, unit_size> MAT<T, unit_size>::operator- () const {
        MAT<T, unit_size> ans = *this;
        ans.neg_();
        return ans;
    }
}

#endif // __QUADTREE_HPP__