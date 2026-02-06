#include "math/sparse_matrix.h"
#include "core/logger.h"
#include <algorithm>
#include <map>

namespace fem {

// ═══ SparseMatrixCOO ═══
void SparseMatrixCOO::print(const std::string& name) const {
    if (!name.empty()) {
        FEM_INFO("SparseMatrixCOO " + name + " (" + std::to_string(rows_) + 
                 "x" + std::to_string(cols_) + ", nnz=" + std::to_string(nnz()) + "):");
    } else {
        FEM_INFO("SparseMatrixCOO (" + std::to_string(rows_) + 
                 "x" + std::to_string(cols_) + ", nnz=" + std::to_string(nnz()) + "):");
    }
    
    std::size_t print_limit = 10;
    for (std::size_t k = 0; k < std::min(nnz(), print_limit); ++k) {
        FEM_INFO("  (" + std::to_string(row_indices_[k]) + ", " + 
                 std::to_string(col_indices_[k]) + ") = " + fmt_sci(values_[k]));
    }
    if (nnz() > print_limit) {
        FEM_INFO("  ... (" + std::to_string(nnz() - print_limit) + " more)");
    }
}

// ═══ SparseMatrixCSR ═══
// ── std::vector 版本 (兼容旧代码) ──
void SparseMatrixCSR::matvec(const Real* x, Real* y) const {
    for (std::size_t i = 0; i < rows_; ++i) {
        Real sum = 0.0;
        for (Index k = row_ptr_[i]; k < row_ptr_[i + 1]; ++k) {
            sum += values_[k] * x[col_indices_[k]];
        }
        y[i] = sum;
    }
}

Vector SparseMatrixCSR::matvec(const Vector& x) const {
    if (cols_ != x.size()) {
        throw std::runtime_error("Matrix-vector size mismatch in CSR matvec");
    }
    
    Vector y(rows_, 0.0);
    for (std::size_t i = 0; i < rows_; ++i) {
        Real sum = 0.0;
        for (Index k = row_ptr_[i]; k < row_ptr_[i + 1]; ++k) {
            sum += values_[k] * x[col_indices_[k]];
        }
        y[i] = sum;
    }
    return y;
}

void SparseMatrixCSR::print(const std::string& name) const {
    if (!name.empty()) {
        FEM_INFO("SparseMatrixCSR " + name + " (" + std::to_string(rows_) + 
                 "x" + std::to_string(cols_) + ", nnz=" + std::to_string(nnz()) + "):");
    } else {
        FEM_INFO("SparseMatrixCSR (" + std::to_string(rows_) + 
                 "x" + std::to_string(cols_) + ", nnz=" + std::to_string(nnz()) + "):");
    }
    
    std::size_t print_rows = std::min(rows_, std::size_t(5));
    for (std::size_t i = 0; i < print_rows; ++i) {
        std::string row_str = "  row " + std::to_string(i) + ": ";
        Index row_begin = row_ptr_[i];
        Index row_end = row_ptr_[i + 1];
        std::size_t count = 0;
        for (Index k = row_begin; k < row_end && count < 5; ++k, ++count) {
            row_str += "(" + std::to_string(col_indices_[k]) + "," + 
                       fmt_sci(values_[k]) + ") ";
        }
        if (row_end - row_begin > 5) {
            row_str += "...";
        }
        FEM_INFO(row_str);
    }
    if (rows_ > print_rows) {
        FEM_INFO("  ... (" + std::to_string(rows_ - print_rows) + " more rows)");
    }
}

// ═══ SparseMatrixCSC ═══
void SparseMatrixCSC::print(const std::string& name) const {
    if (!name.empty()) {
        FEM_INFO("SparseMatrixCSC " + name + " (" + std::to_string(rows_) + 
                 "x" + std::to_string(cols_) + ", nnz=" + std::to_string(nnz()) + "):");
    } else {
        FEM_INFO("SparseMatrixCSC (" + std::to_string(rows_) + 
                 "x" + std::to_string(cols_) + ", nnz=" + std::to_string(nnz()) + "):");
    }
    FEM_INFO("  CSC format: col_ptr.size=" + std::to_string(col_ptr_.size()) +
             ", row_indices.size=" + std::to_string(row_indices_.size()));
}

// ═══ 格式转换 ═══
SparseMatrixCSR coo_to_csr(const SparseMatrixCOO& coo) {
    std::size_t rows = coo.rows();
    std::size_t cols = coo.cols();
    std::size_t nnz = coo.nnz();
    
    if (nnz == 0) {
        SparseMatrixCSR csr(rows, cols);
        std::vector<Index> row_ptr(rows + 1, 0);
        csr.set_data(rows, cols, std::move(row_ptr), {}, {});
        return csr;
    }
    
    // 创建 (row, col, val) 三元组列表并排序
    struct Entry {
        Index row, col;
        Real val;
        bool operator<(const Entry& other) const {
            if (row != other.row) return row < other.row;
            return col < other.col;
        }
    };
    
    std::vector<Entry> entries;
    entries.reserve(nnz);
    for (std::size_t k = 0; k < nnz; ++k) {
        entries.push_back({coo.row_indices()[k], coo.col_indices()[k], coo.values()[k]});
    }
    std::sort(entries.begin(), entries.end());
    
    // 合并重复项 (同一位置的值累加)
    std::vector<Entry> merged;
    merged.push_back(entries[0]);
    for (std::size_t k = 1; k < entries.size(); ++k) {
        if (entries[k].row == merged.back().row && 
            entries[k].col == merged.back().col) {
            merged.back().val += entries[k].val;
        } else {
            merged.push_back(entries[k]);
        }
    }
    
    // 构建 CSR
    std::vector<Index> row_ptr(rows + 1, 0);
    std::vector<Index> col_indices;
    std::vector<Real> values;
    
    col_indices.reserve(merged.size());
    values.reserve(merged.size());
    
    Index current_row = 0;
    for (const auto& e : merged) {
        // 填充空行
        while (current_row < e.row) {
            row_ptr[current_row + 1] = static_cast<Index>(col_indices.size());
            ++current_row;
        }
        
        col_indices.push_back(e.col);
        values.push_back(e.val);
    }
    
    // 填充剩余空行
    while (current_row < rows) {
        row_ptr[current_row + 1] = static_cast<Index>(col_indices.size());
        ++current_row;
    }
    
    SparseMatrixCSR csr(rows, cols);
    csr.set_data(rows, cols, std::move(row_ptr), std::move(col_indices), std::move(values));
    return csr;
}

SparseMatrixCSC coo_to_csc(const SparseMatrixCOO& coo) {
    std::size_t rows = coo.rows();
    std::size_t cols = coo.cols();
    std::size_t nnz = coo.nnz();
    
    if (nnz == 0) {
        SparseMatrixCSC csc(rows, cols);
        std::vector<Index> col_ptr(cols + 1, 0);
        csc.set_data(rows, cols, std::move(col_ptr), {}, {});
        return csc;
    }
    
    // 创建 (col, row, val) 三元组列表并排序
    struct Entry {
        Index col, row;
        Real val;
        bool operator<(const Entry& other) const {
            if (col != other.col) return col < other.col;
            return row < other.row;
        }
    };
    
    std::vector<Entry> entries;
    entries.reserve(nnz);
    for (std::size_t k = 0; k < nnz; ++k) {
        entries.push_back({coo.col_indices()[k], coo.row_indices()[k], coo.values()[k]});
    }
    std::sort(entries.begin(), entries.end());
    
    // 合并重复项
    std::vector<Entry> merged;
    merged.push_back(entries[0]);
    for (std::size_t k = 1; k < entries.size(); ++k) {
        if (entries[k].col == merged.back().col && 
            entries[k].row == merged.back().row) {
            merged.back().val += entries[k].val;
        } else {
            merged.push_back(entries[k]);
        }
    }
    
    // 构建 CSC
    std::vector<Index> col_ptr(cols + 1, 0);
    std::vector<Index> row_indices;
    std::vector<Real> values;
    
    row_indices.reserve(merged.size());
    values.reserve(merged.size());
    
    Index current_col = 0;
    for (const auto& e : merged) {
        while (current_col < e.col) {
            col_ptr[current_col + 1] = static_cast<Index>(row_indices.size());
            ++current_col;
        }
        
        row_indices.push_back(e.row);
        values.push_back(e.val);
    }
    
    while (current_col < cols) {
        col_ptr[current_col + 1] = static_cast<Index>(row_indices.size());
        ++current_col;
    }
    
    SparseMatrixCSC csc(rows, cols);
    csc.set_data(rows, cols, std::move(col_ptr), std::move(row_indices), std::move(values));
    return csc;
}

SparseMatrixCSR csc_to_csr(const SparseMatrixCSC& csc) {
    // 先转 COO 再转 CSR (简化实现)
    SparseMatrixCOO coo(csc.rows(), csc.cols());
    const auto& col_ptr = csc.col_ptr();
    const auto& row_indices = csc.row_indices();
    const auto& values = csc.values();
    
    for (std::size_t j = 0; j < csc.cols(); ++j) {
        for (Index k = col_ptr[j]; k < col_ptr[j + 1]; ++k) {
            coo.add(row_indices[k], j, values[k]);
        }
    }
    
    return coo_to_csr(coo);
}

SparseMatrixCSC csr_to_csc(const SparseMatrixCSR& csr) {
    // 先转 COO 再转 CSC
    SparseMatrixCOO coo(csr.rows(), csr.cols());
    const auto& row_ptr = csr.row_ptr();
    const auto& col_indices = csr.col_indices();
    const auto& values = csr.values();
    
    for (std::size_t i = 0; i < csr.rows(); ++i) {
        for (Index k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
            coo.add(i, col_indices[k], values[k]);
        }
    }
    
    return coo_to_csc(coo);
}

// ── CSR → COO 转换 ──
SparseMatrixCOO csr_to_coo(const SparseMatrixCSR& csr) {
    SparseMatrixCOO coo(csr.rows(), csr.cols());
    
    for (std::size_t i = 0; i < csr.rows(); ++i) {
        Index start = csr.row_ptr()[i];
        Index end = csr.row_ptr()[i + 1];
        
        for (Index k = start; k < end; ++k) {
            Index j = csr.col_indices()[k];
            Real val = csr.values()[k];
            coo.add(i, j, val);
        }
    }
    
    return coo;
}

}  // namespace fem
