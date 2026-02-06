#pragma once

#include "core/types.h"
#include <vector>

namespace fem {

// ── COO 格式 (装配阶段: 随机 += 快) ──
struct COOMatrix {
    std::size_t rows{0};
    std::size_t cols{0};
    std::vector<Index> row_idx;
    std::vector<Index> col_idx;
    std::vector<Real>  values;

    void add(Index i, Index j, Real val) {
        row_idx.push_back(i);
        col_idx.push_back(j);
        values.push_back(val);
    }
};

// ── CSR 格式 (求解阶段: 按行遍历快) ──
struct CSRMatrix {
    std::size_t rows{0};
    std::vector<Index> row_ptr;   // 长度 rows+1
    std::vector<Index> col_idx;
    std::vector<Real>  values;

    // y = A * x
    void matvec(const Real* x, Real* y) const;
};

// COO → CSR 转换 (重复项自动求和)
[[nodiscard]] CSRMatrix coo_to_csr(const COOMatrix& coo);

}  // namespace fem
