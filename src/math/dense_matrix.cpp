#include "math/dense_matrix.h"
#include "core/logger.h"
#include <cmath>
#include <algorithm>

namespace fem {

void DenseMatrix::identity() {
    if (!is_square()) {
        throw std::runtime_error("Identity matrix must be square");
    }
    zero();
    for (std::size_t i = 0; i < rows_; ++i) {
        (*this)(i, i) = 1.0;
    }
}

Vector DenseMatrix::get_row(std::size_t i) const {
    if (i >= rows_) {
        throw std::out_of_range("Row index out of range");
    }
    Vector v(cols_);
    for (std::size_t j = 0; j < cols_; ++j) {
        v[j] = (*this)(i, j);
    }
    return v;
}

Vector DenseMatrix::get_col(std::size_t j) const {
    if (j >= cols_) {
        throw std::out_of_range("Column index out of range");
    }
    Vector v(rows_);
    for (std::size_t i = 0; i < rows_; ++i) {
        v[i] = (*this)(i, j);
    }
    return v;
}

void DenseMatrix::set_row(std::size_t i, const Vector& v) {
    if (i >= rows_ || v.size() != cols_) {
        throw std::runtime_error("Row size mismatch");
    }
    for (std::size_t j = 0; j < cols_; ++j) {
        (*this)(i, j) = v[j];
    }
}

void DenseMatrix::set_col(std::size_t j, const Vector& v) {
    if (j >= cols_ || v.size() != rows_) {
        throw std::runtime_error("Column size mismatch");
    }
    for (std::size_t i = 0; i < rows_; ++i) {
        (*this)(i, j) = v[i];
    }
}

DenseMatrix DenseMatrix::operator+(const DenseMatrix& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        throw std::runtime_error("Matrix size mismatch in addition");
    }
    DenseMatrix result(rows_, cols_);
    for (std::size_t i = 0; i < data_.size(); ++i) {
        result.data_[i] = data_[i] + other.data_[i];
    }
    return result;
}

DenseMatrix DenseMatrix::operator-(const DenseMatrix& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        throw std::runtime_error("Matrix size mismatch in subtraction");
    }
    DenseMatrix result(rows_, cols_);
    for (std::size_t i = 0; i < data_.size(); ++i) {
        result.data_[i] = data_[i] - other.data_[i];
    }
    return result;
}

DenseMatrix DenseMatrix::operator*(Real scalar) const {
    DenseMatrix result(rows_, cols_);
    for (std::size_t i = 0; i < data_.size(); ++i) {
        result.data_[i] = data_[i] * scalar;
    }
    return result;
}

DenseMatrix DenseMatrix::operator*(const DenseMatrix& other) const {
    if (cols_ != other.rows_) {
        throw std::runtime_error("Matrix size mismatch in multiplication");
    }
    DenseMatrix result(rows_, other.cols_, 0.0);
    for (std::size_t i = 0; i < rows_; ++i) {
        for (std::size_t j = 0; j < other.cols_; ++j) {
            Real sum = 0.0;
            for (std::size_t k = 0; k < cols_; ++k) {
                sum += (*this)(i, k) * other(k, j);
            }
            result(i, j) = sum;
        }
    }
    return result;
}

DenseMatrix& DenseMatrix::operator+=(const DenseMatrix& other) {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        throw std::runtime_error("Matrix size mismatch in addition");
    }
    for (std::size_t i = 0; i < data_.size(); ++i) {
        data_[i] += other.data_[i];
    }
    return *this;
}

DenseMatrix& DenseMatrix::operator-=(const DenseMatrix& other) {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        throw std::runtime_error("Matrix size mismatch in subtraction");
    }
    for (std::size_t i = 0; i < data_.size(); ++i) {
        data_[i] -= other.data_[i];
    }
    return *this;
}

DenseMatrix& DenseMatrix::operator*=(Real scalar) {
    for (Real& x : data_) {
        x *= scalar;
    }
    return *this;
}

Vector DenseMatrix::matvec(const Vector& v) const {
    if (cols_ != v.size()) {
        throw std::runtime_error("Matrix-vector size mismatch");
    }
    Vector result(rows_, 0.0);
    for (std::size_t i = 0; i < rows_; ++i) {
        Real sum = 0.0;
        for (std::size_t j = 0; j < cols_; ++j) {
            sum += (*this)(i, j) * v[j];
        }
        result[i] = sum;
    }
    return result;
}

DenseMatrix DenseMatrix::transpose() const {
    DenseMatrix result(cols_, rows_);
    for (std::size_t i = 0; i < rows_; ++i) {
        for (std::size_t j = 0; j < cols_; ++j) {
            result(j, i) = (*this)(i, j);
        }
    }
    return result;
}

void DenseMatrix::transpose_inplace() {
    if (!is_square()) {
        throw std::runtime_error("In-place transpose requires square matrix");
    }
    for (std::size_t i = 0; i < rows_; ++i) {
        for (std::size_t j = i + 1; j < cols_; ++j) {
            std::swap((*this)(i, j), (*this)(j, i));
        }
    }
}

bool DenseMatrix::is_symmetric(Real tol) const {
    if (!is_square()) return false;
    for (std::size_t i = 0; i < rows_; ++i) {
        for (std::size_t j = i + 1; j < cols_; ++j) {
            if (std::abs((*this)(i, j) - (*this)(j, i)) > tol) {
                return false;
            }
        }
    }
    return true;
}

Real DenseMatrix::norm_frobenius() const {
    Real sum = 0.0;
    for (Real x : data_) {
        sum += x * x;
    }
    return std::sqrt(sum);
}

Real DenseMatrix::norm_inf() const {
    Real max_row_sum = 0.0;
    for (std::size_t i = 0; i < rows_; ++i) {
        Real row_sum = 0.0;
        for (std::size_t j = 0; j < cols_; ++j) {
            row_sum += std::abs((*this)(i, j));
        }
        max_row_sum = std::max(max_row_sum, row_sum);
    }
    return max_row_sum;
}

void DenseMatrix::print(const std::string& name) const {
    if (!name.empty()) {
        FEM_INFO("DenseMatrix " + name + " (" + std::to_string(rows_) + 
                 "x" + std::to_string(cols_) + "):");
    } else {
        FEM_INFO("DenseMatrix (" + std::to_string(rows_) + 
                 "x" + std::to_string(cols_) + "):");
    }
    
    std::size_t print_rows = std::min(rows_, std::size_t(5));
    std::size_t print_cols = std::min(cols_, std::size_t(5));
    
    for (std::size_t i = 0; i < print_rows; ++i) {
        std::string row_str = "  [" + std::to_string(i) + "] ";
        for (std::size_t j = 0; j < print_cols; ++j) {
            row_str += fmt_sci((*this)(i, j)) + " ";
        }
        if (cols_ > print_cols) {
            row_str += "...";
        }
        FEM_INFO(row_str);
    }
    if (rows_ > print_rows) {
        FEM_INFO("  ... (" + std::to_string(rows_ - print_rows) + " more rows)");
    }
}

}  // namespace fem
