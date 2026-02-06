#pragma once

#include "math/vector.h"
#include "core/types.h"
#include <vector>
#include <tuple>

namespace fem {

/**
 * SparseMatrixCOO: COO (Coordinate) 格式稀疏矩阵
 * 
 * 存储 (row, col, value) 三元组
 * 适合矩阵装配阶段
 */
class SparseMatrixCOO {
public:
    // ═══ 构造 ═══
    SparseMatrixCOO() : rows_(0), cols_(0) {}
    SparseMatrixCOO(std::size_t rows, std::size_t cols)
        : rows_(rows), cols_(cols) {}
    
    // ═══ 大小 ═══
    std::size_t rows() const { return rows_; }
    std::size_t cols() const { return cols_; }
    std::size_t nnz() const { return row_indices_.size(); }  // 非零元数量
    
    void set_size(std::size_t rows, std::size_t cols) {
        rows_ = rows;
        cols_ = cols;
    }
    
    // ═══ 添加元素 ═══
    void add(Index i, Index j, Real value) {
        row_indices_.push_back(i);
        col_indices_.push_back(j);
        values_.push_back(value);
    }
    
    void reserve(std::size_t n) {
        row_indices_.reserve(n);
        col_indices_.reserve(n);
        values_.reserve(n);
    }
    
    void clear() {
        row_indices_.clear();
        col_indices_.clear();
        values_.clear();
    }
    
    // ═══ 访问数据 ═══
    const std::vector<Index>& row_indices() const { return row_indices_; }
    const std::vector<Index>& col_indices() const { return col_indices_; }
    const std::vector<Real>& values() const { return values_; }
    
    // ═══ 工具 ═══
    void print(const std::string& name = "") const;
    
private:
    std::size_t rows_;
    std::size_t cols_;
    std::vector<Index> row_indices_;
    std::vector<Index> col_indices_;
    std::vector<Real> values_;
};

/**
 * SparseMatrixCSR: CSR (Compressed Sparse Row) 格式
 * 
 * 行压缩存储，适合矩阵-向量乘法
 * row_ptr[i] ~ row_ptr[i+1]: 第i行的非零元在 col_indices 和 values 中的范围
 */
class SparseMatrixCSR {
public:
    // ═══ 构造 ═══
    SparseMatrixCSR() : rows_(0), cols_(0) {}
    SparseMatrixCSR(std::size_t rows, std::size_t cols)
        : rows_(rows), cols_(cols) {}
    
    // ═══ 大小 ═══
    std::size_t rows() const { return rows_; }
    std::size_t cols() const { return cols_; }
    std::size_t nnz() const { return values_.size(); }
    
    // ═══ 设置数据 ═══
    void set_data(std::size_t rows, std::size_t cols,
                  std::vector<Index>&& row_ptr,
                  std::vector<Index>&& col_indices,
                  std::vector<Real>&& values) {
        rows_ = rows;
        cols_ = cols;
        row_ptr_ = std::move(row_ptr);
        col_indices_ = std::move(col_indices);
        values_ = std::move(values);
    }
    
    // ═══ 访问数据 ═══
    const std::vector<Index>& row_ptr() const { return row_ptr_; }
    const std::vector<Index>& col_indices() const { return col_indices_; }
    const std::vector<Real>& values() const { return values_; }
    
    std::vector<Index>& row_ptr() { return row_ptr_; }
    std::vector<Index>& col_indices() { return col_indices_; }
    std::vector<Real>& values() { return values_; }
    
    // ═══ 矩阵-向量乘法 ═══
    Vector matvec(const Vector& x) const;
    Vector operator*(const Vector& x) const { return matvec(x); }
    
    // ═══ 工具 ═══
    void print(const std::string& name = "") const;
    
private:
    std::size_t rows_;
    std::size_t cols_;
    std::vector<Index> row_ptr_;      // 长度 rows+1
    std::vector<Index> col_indices_;  // 长度 nnz
    std::vector<Real> values_;        // 长度 nnz
};

/**
 * SparseMatrixCSC: CSC (Compressed Sparse Column) 格式
 * 
 * 列压缩存储
 */
class SparseMatrixCSC {
public:
    // ═══ 构造 ═══
    SparseMatrixCSC() : rows_(0), cols_(0) {}
    SparseMatrixCSC(std::size_t rows, std::size_t cols)
        : rows_(rows), cols_(cols) {}
    
    // ═══ 大小 ═══
    std::size_t rows() const { return rows_; }
    std::size_t cols() const { return cols_; }
    std::size_t nnz() const { return values_.size(); }
    
    // ═══ 设置数据 ═══
    void set_data(std::size_t rows, std::size_t cols,
                  std::vector<Index>&& col_ptr,
                  std::vector<Index>&& row_indices,
                  std::vector<Real>&& values) {
        rows_ = rows;
        cols_ = cols;
        col_ptr_ = std::move(col_ptr);
        row_indices_ = std::move(row_indices);
        values_ = std::move(values);
    }
    
    // ═══ 访问数据 ═══
    const std::vector<Index>& col_ptr() const { return col_ptr_; }
    const std::vector<Index>& row_indices() const { return row_indices_; }
    const std::vector<Real>& values() const { return values_; }
    
    // ═══ 工具 ═══
    void print(const std::string& name = "") const;
    
private:
    std::size_t rows_;
    std::size_t cols_;
    std::vector<Index> col_ptr_;      // 长度 cols+1
    std::vector<Index> row_indices_;  // 长度 nnz
    std::vector<Real> values_;        // 长度 nnz
};

// ═══ 格式转换 ═══
SparseMatrixCSR coo_to_csr(const SparseMatrixCOO& coo);
SparseMatrixCSC coo_to_csc(const SparseMatrixCOO& coo);
SparseMatrixCSR csc_to_csr(const SparseMatrixCSC& csc);
SparseMatrixCSC csr_to_csc(const SparseMatrixCSR& csr);
SparseMatrixCOO csr_to_coo(const SparseMatrixCSR& csr);

}  // namespace fem
