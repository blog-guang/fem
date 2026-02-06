#include "assembly/boundary_condition.h"
#include <algorithm>
#include <unordered_set>

namespace fem {

void apply_dirichlet(CSRMatrix&    K,
                     Vector&       F,
                     const Mesh&   mesh,
                     const BoundaryCondition& bc,
                     std::size_t   dofs_per_node)
{
    if (bc.type != BCType::Dirichlet) return;

    const auto& bc_nodes = mesh.boundary_nodes(bc.boundary_name);
    
    // 计算需要固定的全局自由度
    std::vector<Index> fixed_dofs;
    for (Index node : bc_nodes) {
        if (bc.component == -1) {
            // 固定所有分量
            for (std::size_t dof = 0; dof < dofs_per_node; ++dof) {
                fixed_dofs.push_back(node * dofs_per_node + dof);
            }
        } else if (bc.component >= 0 && bc.component < static_cast<int>(dofs_per_node)) {
            // 只固定特定分量
            fixed_dofs.push_back(node * dofs_per_node + static_cast<Index>(bc.component));
        }
    }

    std::unordered_set<Index> fixed_dof_set(fixed_dofs.begin(), fixed_dofs.end());

    // 1. 对固定自由度: 清除行, 对角线置 1, F = value
    for (Index dof : fixed_dofs) {
        if (dof >= K.rows) continue;  // 边界检查
        
        for (Index j = K.row_ptr[dof]; j < K.row_ptr[dof + 1]; ++j) {
            K.values[j] = (K.col_idx[j] == dof) ? 1.0 : 0.0;
        }
        F[dof] = bc.value;
    }

    // 2. 对非固定自由度: 清除涉及固定自由度的列, 将贡献移到 F
    for (std::size_t i = 0; i < K.rows; ++i) {
        if (fixed_dof_set.count(i)) continue;  // 已处理

        for (Index j = K.row_ptr[i]; j < K.row_ptr[i + 1]; ++j) {
            if (fixed_dof_set.count(K.col_idx[j])) {
                F[i] -= K.values[j] * bc.value;
                K.values[j] = 0.0;
            }
        }
    }
}

void apply_boundary_conditions(CSRMatrix&                    K,
                               Vector&                       F,
                               const Mesh&                   mesh,
                               const std::vector<BoundaryCondition>& bcs,
                               std::size_t                   dofs_per_node)
{
    for (const auto& bc : bcs) {
        apply_dirichlet(K, F, mesh, bc, dofs_per_node);
    }
}

}  // namespace fem
