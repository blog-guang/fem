#include "assembly/assembler.h"

namespace fem {

Assembler::Assembler(const Mesh& mesh, const ElementBase& elem, std::size_t dofs_per_node)
    : mesh_(mesh), elem_(elem), dofs_per_node_(dofs_per_node) {}

void Assembler::assemble(ElementStiffnessFn ke_fn,
                          ElementLoadFn      fe_fn,
                          COOMatrix&         K,
                          Vector&            F,
                          void*              ctx) const
{
    std::size_t num_nodes = mesh_.num_nodes();
    std::size_t total_dofs = num_nodes * dofs_per_node_;
    
    K.rows = total_dofs;
    K.cols = total_dofs;
    F.assign(total_dofs, 0.0);

    // 临时缓冲区 (最大支持 MAX_NODES * dofs_per_node 自由度)
    constexpr std::size_t MAX_LOCAL_SIZE = MAX_NODES * 8; // assuming max dofs_per_node = 8
    Real Ke[MAX_LOCAL_SIZE * MAX_LOCAL_SIZE];
    Real Fe[MAX_LOCAL_SIZE];
    Vec3 coords[MAX_NODES];

    for (std::size_t c = 0; c < mesh_.num_cells(); ++c) {
        const Cell& cell = mesh_.cell(c);
        std::size_t n    = cell.num_nodes;

        // 提取单元节点坐标
        for (std::size_t i = 0; i < n; ++i) {
            coords[i] = mesh_.coords(cell.node(i));
        }

        // 计算局部自由度总数
        std::size_t local_size = n * dofs_per_node_;

        // 单元刚度矩阵
        if (ke_fn) {
            ke_fn(coords, n, dofs_per_node_, Ke, ctx);
            // 装配到全局 K
            for (std::size_t i = 0; i < n; ++i) {
                for (std::size_t j = 0; j < n; ++j) {
                    for (std::size_t di = 0; di < dofs_per_node_; ++di) {
                        for (std::size_t dj = 0; dj < dofs_per_node_; ++dj) {
                            Index gi = cell.node(i) * dofs_per_node_ + di;  // 全局 DOF i
                            Index gj = cell.node(j) * dofs_per_node_ + dj;  // 全局 DOF j
                            K.add(gi, gj, Ke[(i * dofs_per_node_ + di) * local_size + (j * dofs_per_node_ + dj)]);
                        }
                    }
                }
            }
        }

        // 单元载荷向量
        if (fe_fn) {
            fe_fn(coords, n, dofs_per_node_, Fe, ctx);
            // 装配到全局 F
            for (std::size_t i = 0; i < n; ++i) {
                for (std::size_t di = 0; di < dofs_per_node_; ++di) {
                    Index gi = cell.node(i) * dofs_per_node_ + di;  // 全局 DOF
                    F[gi] += Fe[i * dofs_per_node_ + di];
                }
            }
        }
    }
}

}  // namespace fem
