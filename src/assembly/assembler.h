#pragma once

#include "core/types.h"
#include "mesh/mesh.h"
#include "element/element_base.h"
#include "assembly/sparse_matrix.h"

namespace fem {

// ── 装配器 ──
// 遍历网格所有单元, 调用用户提供的单元矩阵/向量函数, 装配到全局 COO + F
//
// 用法:
//   Assembler asm(mesh, element, dofs_per_node);
//   asm.assemble(element_stiffness_fn, element_load_fn, K, F);
//

// 单元刚度矩阵回调: 给定单元节点坐标、节点数、自由度数 → 填入 Ke[(n*dof)x(n*dof)]
using ElementStiffnessFn =
    void(*)(const Vec3* coords, std::size_t n_nodes, std::size_t dofs_per_node, Real* Ke, void* ctx);

// 单元载荷向量回调: 给定单元节点坐标、节点数、自由度数 → 填入 Fe[n*dof]
using ElementLoadFn =
    void(*)(const Vec3* coords, std::size_t n_nodes, std::size_t dofs_per_node, Real* Fe, void* ctx);

class Assembler {
public:
    explicit Assembler(const Mesh& mesh, const ElementBase& elem, std::size_t dofs_per_node = 1);

    // 装配全局刚度矩阵 K (COO) + 载荷向量 F
    void assemble(ElementStiffnessFn ke_fn,
                  ElementLoadFn      fe_fn,
                  COOMatrix&         K,
                  Vector&            F,
                  void*              ctx = nullptr) const;

private:
    const Mesh&        mesh_;
    const ElementBase& elem_;
    std::size_t        dofs_per_node_;
};

}  // namespace fem
