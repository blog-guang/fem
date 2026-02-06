#include "assembly/sparse_matrix.h"
#include <algorithm>
#include <cstring>

namespace fem {

void CSRMatrix::matvec(const Real* x, Real* y) const {
    std::memset(y, 0, rows * sizeof(Real));
    for (std::size_t i = 0; i < rows; ++i) {
        Real sum = 0.0;
        for (Index j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
            sum += values[j] * x[col_idx[j]];
        }
        y[i] = sum;
    }
}

CSRMatrix coo_to_csr(const COOMatrix& coo) {
    std::size_t nnz = coo.values.size();

    // 1. 按 (row, col) 排序
    std::vector<std::size_t> order(nnz);
    for (std::size_t i = 0; i < nnz; ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        if (coo.row_idx[a] != coo.row_idx[b])
            return coo.row_idx[a] < coo.row_idx[b];
        return coo.col_idx[a] < coo.col_idx[b];
    });

    // 2. 遍历排序后数据, 同一 (row, col) 求和 → merged arrays
    std::vector<Index> m_row, m_col;
    std::vector<Real>  m_val;
    m_row.reserve(nnz);
    m_col.reserve(nnz);
    m_val.reserve(nnz);

    for (std::size_t k = 0; k < nnz; ++k) {
        Index r = coo.row_idx[order[k]];
        Index c = coo.col_idx[order[k]];
        Real  v = coo.values[order[k]];

        if (!m_row.empty() && r == m_row.back() && c == m_col.back()) {
            m_val.back() += v;   // 重复项累加
        } else {
            m_row.push_back(r);
            m_col.push_back(c);
            m_val.push_back(v);
        }
    }

    // 3. 构建 row_ptr (从 m_row 统计每行元素数)
    CSRMatrix csr;
    csr.rows    = coo.rows;
    csr.row_ptr.assign(coo.rows + 1, 0);
    for (Index r : m_row) {
        csr.row_ptr[r + 1]++;
    }
    for (std::size_t i = 1; i <= coo.rows; ++i) {
        csr.row_ptr[i] += csr.row_ptr[i - 1];
    }

    csr.col_idx = std::move(m_col);
    csr.values  = std::move(m_val);
    return csr;
}

}  // namespace fem
