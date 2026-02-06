#pragma once

#include "core/types.h"
#include "mesh/mesh.h"
#include "assembly/sparse_matrix.h"

namespace fem {

// ── 边界条件 ──
struct BoundaryCondition {
    BCType      type;             // Dirichlet / Neumann
    std::string boundary_name;   // 对应 Mesh 中的边界标记
    Real        value{0.0};      // 常数值
    int         component{-1};   // 分量索引 (-1=all components for scalar, 0=x, 1=y, etc.)
};

// ── 施加 Dirichlet BC 到 CSR 矩阵 + 载荷向量 ──
// 对 BC 节点: 将该行置为单位行 (对角线=1, 其他=0), F[i] = value
// 同时修改其他行中涉及该节点的列 (保持对称性)
void apply_dirichlet(CSRMatrix&    K,
                     Vector&       F,
                     const Mesh&   mesh,
                     const BoundaryCondition& bc,
                     std::size_t   dofs_per_node = 1);

// 批量施加
void apply_boundary_conditions(CSRMatrix&                    K,
                               Vector&                       F,
                               const Mesh&                   mesh,
                               const std::vector<BoundaryCondition>& bcs,
                               std::size_t                   dofs_per_node = 1);

}  // namespace fem
