#pragma once

#include "math/vector.h"
#include "core/types.h"
#include <vector>

namespace fem {

/**
 * DenseMatrix: 密集矩阵类
 * 
 * 行优先存储 (row-major order)
 * data[i*cols + j] = A(i,j)
 */
class DenseMatrix {
public:
    // ═══ 构造 ═══
    DenseMatrix() : rows_(0), cols_(0) {}
    DenseMatrix(std::size_t rows, std::size_t cols, Real value = 0.0)
        : rows_(rows), cols_(cols), data_(rows * cols, value) {}
    
    // ═══ 大小与访问 ═══
    std::size_t rows() const { return rows_; }
    std::size_t cols() const { return cols_; }
    std::size_t size() const { return rows_ * cols_; }
    
    void resize(std::size_t rows, std::size_t cols, Real value = 0.0) {
        rows_ = rows;
        cols_ = cols;
        data_.resize(rows * cols, value);
    }
    
    void clear() {
        rows_ = 0;
        cols_ = 0;
        data_.clear();
    }
    
    // 元素访问: A(i, j)
    Real& operator()(std::size_t i, std::size_t j) {
        return data_[i * cols_ + j];
    }
    const Real& operator()(std::size_t i, std::size_t j) const {
        return data_[i * cols_ + j];
    }
    
    Real& at(std::size_t i, std::size_t j) {
        if (i >= rows_ || j >= cols_) {
            throw std::out_of_range("Matrix index out of range");
        }
        return data_[i * cols_ + j];
    }
    const Real& at(std::size_t i, std::size_t j) const {
        if (i >= rows_ || j >= cols_) {
            throw std::out_of_range("Matrix index out of range");
        }
        return data_[i * cols_ + j];
    }
    
    Real* data() { return data_.data(); }
    const Real* data() const { return data_.data(); }
    
    // ═══ 赋值 ═══
    void fill(Real value) { std::fill(data_.begin(), data_.end(), value); }
    void zero() { fill(0.0); }
    void identity();  // 设为单位矩阵
    
    // ═══ 行列操作 ═══
    Vector get_row(std::size_t i) const;
    Vector get_col(std::size_t j) const;
    void set_row(std::size_t i, const Vector& v);
    void set_col(std::size_t j, const Vector& v);
    
    // ═══ 矩阵运算 ═══
    // A + B
    DenseMatrix operator+(const DenseMatrix& other) const;
    // A - B
    DenseMatrix operator-(const DenseMatrix& other) const;
    // a * A
    DenseMatrix operator*(Real scalar) const;
    // A * B (矩阵乘法)
    DenseMatrix operator*(const DenseMatrix& other) const;
    
    // A += B
    DenseMatrix& operator+=(const DenseMatrix& other);
    // A -= B
    DenseMatrix& operator-=(const DenseMatrix& other);
    // A *= a
    DenseMatrix& operator*=(Real scalar);
    
    // A * v (矩阵-向量乘法)
    Vector matvec(const Vector& v) const;
    Vector operator*(const Vector& v) const { return matvec(v); }
    
    // 转置
    DenseMatrix transpose() const;
    void transpose_inplace();  // 仅方阵
    
    // ═══ 性质 ═══
    bool is_square() const { return rows_ == cols_; }
    bool is_symmetric(Real tol = 1e-12) const;
    
    // 范数
    Real norm_frobenius() const;  // Frobenius 范数
    Real norm_inf() const;        // 无穷范数 (最大行和)
    
    // ═══ 工具函数 ═══
    void print(const std::string& name = "") const;
    
private:
    std::size_t rows_;
    std::size_t cols_;
    std::vector<Real> data_;  // 行优先存储
};

// ═══ 全局运算符 ═══
inline DenseMatrix operator*(Real scalar, const DenseMatrix& A) {
    return A * scalar;
}

}  // namespace fem
