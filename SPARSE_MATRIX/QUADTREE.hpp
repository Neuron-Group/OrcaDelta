
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
        const TREENODE<T> & operator = (const TREENODE<T> &b);
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
    const TREENODE<T> & TREENODE<T>::operator = (const TREENODE<T> &b) {

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

        return b;
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

    // 矩阵类
    template<class T, size_t unit_size>
    class MAT {
        size_t rows, cols;

        TREENODE<Eigen::Matrix<T, unit_size, unit_size>> * ROOT;

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

        public:
        MAT(size_t rows = 1, size_t cols = 1);
        const MAT<T, unit_size> & operator = (const MAT<T, unit_size> &b);
        MAT(const MAT &b);
        ~MAT();

        bool isBlank () const {
            return (this->ROOT == nullptr) ? true : false;
        }

        bool insert(T data, size_t i, size_t j){
            check_idx_with_error(i, j);
            if (data == static_cast<T>(0)) return true;
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

        T getVal(size_t i, size_t j) {
            check_idx_with_error(i, j);
            // 索引到块
            TREENODE<Eigen::Matrix<T, unit_size, unit_size>> * ptr = this->ROOT->get(IDX(i), IDX(j));
            
            // 在块中索引出元素
            if (ptr == nullptr) return static_cast<T>(0);
            if (ptr->dataNode == nullptr) return static_cast<T>(0);
            return (*(ptr->dataNode))(IDX_(i), IDX_(j));
        }

        void test() {
            ROOT->test_print(ROOT);
        }

        // 转置
        void transpose_();
        MAT<T, unit_size> transpose() const;
    };

    template<class T, size_t unit_size>
    MAT<T, unit_size>::MAT(size_t rows, size_t cols) {
        this->rows = rows;
        this->cols = cols;
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
    const MAT<T, unit_size> & MAT<T, unit_size>::operator = (const MAT<T, unit_size> &b) {
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
        
        q_ptr_this.push(this->ROOT);
        q_ptr_b.push(b.ROOT);
        while(!q_ptr_this.empty()) {
            TREENODE<Eigen::Matrix<T, unit_size, unit_size>> * ptr_this = q_ptr_this.front();
            TREENODE<Eigen::Matrix<T, unit_size, unit_size>> * ptr_b = q_ptr_b.front();

            for (short i = 0; i < 4; i++){
                if (ptr_b->CHILD[i] == nullptr) {
                    del(ptr_this->CHILD[i]);
                    ptr_this->CHILD[i] = nullptr;
                } else {
                    if (ptr_this->CHILD[i] == nullptr) {
                        ptr_this->CHILD[i] = new TREENODE<Eigen::Matrix<T, unit_size, unit_size>>(*(ptr_b->CHILD[i]));
                    } else {
                        *(ptr_this->CHILD[i]) = *(ptr_b->CHILD[i]);
                    }

                    q_ptr_this.push(ptr_this->CHILD[i]);
                    q_ptr_b.push(ptr_b->CHILD[i]);
                }
            }

            q_ptr_this.pop();
            q_ptr_b.pop();
        }

        return b;
    }

    template<class T, size_t unit_size>
    MAT<T, unit_size>::MAT(const MAT &b) {
        this->rows = b.rows;
        this->cols = b.cols;
        this->ROOT = nullptr;
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
        
        q_ptr_this.push(this->ROOT);
        q_ptr_b.push(b.ROOT);
        while(!q_ptr_this.empty()) {
            TREENODE<Eigen::Matrix<T, unit_size, unit_size>> * ptr_this = q_ptr_this.front();
            TREENODE<Eigen::Matrix<T, unit_size, unit_size>> * ptr_b = q_ptr_b.front();

            for (short i = 0; i < 4; i++){
                if (ptr_b->CHILD[i] == nullptr) {
                    del(ptr_this->CHILD[i]);
                    ptr_this->CHILD[i] = nullptr;
                } else {
                    if (ptr_this->CHILD[i] == nullptr) {
                        ptr_this->CHILD[i] = new TREENODE<Eigen::Matrix<T, unit_size, unit_size>>(*(ptr_b->CHILD[i]));
                    } else {
                        *(ptr_this->CHILD[i]) = *(ptr_b->CHILD[i]);
                    }

                    q_ptr_this.push(ptr_this->CHILD[i]);
                    q_ptr_b.push(ptr_b->CHILD[i]);
                }
            }

            q_ptr_this.pop();
            q_ptr_b.pop();
        }
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
}

#endif // __QUADTREE_HPP__
