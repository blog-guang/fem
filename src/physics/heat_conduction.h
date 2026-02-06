#pragma once

#include "core/types.h"

namespace fem {

// ── 热传导材料属性 ──
struct HeatMaterial {
    Real conductivity{1.0};  // 导热系数 k
    Real source{0.0};        // 热源 Q
};

// ── 热传导单元刚度矩阵 ──
// Ke_ij = k * area * (dNi_dx*dNj_dx + dNi_dy*dNj_dy)
void heat_stiffness(const Vec3* coords, std::size_t n_nodes, std::size_t dofs_per_node, Real* Ke, void* ctx);

// ── 热传导单元载荷向量 ──
// Fe_i = Q * area / 3
void heat_load(const Vec3* coords, std::size_t n_nodes, std::size_t dofs_per_node, Real* Fe, void* ctx);

}  // namespace fem