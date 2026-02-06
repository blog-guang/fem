#pragma once

#include "core/types.h"

namespace fem {

// ── 弹性力学材料属性 (平面应力) ──
struct ElasticMaterial {
    Real E{1.0};           // 弹性模量
    Real poisson{0.0};     // 泊松比
    Real thickness{1.0};   // 厚度 (平面应力)
};

// ── 弹性体力载荷 ──
struct ElasticLoad {
    Real fx{0.0};          // x方向体力
    Real fy{0.0};          // y方向体力
};

// ── 弹性力学上下文 ──
struct ElasticCtx {
    ElasticMaterial mat;
    ElasticLoad load;
};

// ── 弹性力学单元刚度矩阵 ──
// Ke = ∫ B^T * D * B * t * dA
void elasticity_stiffness(const Vec3* coords, std::size_t n_nodes, std::size_t dofs_per_node, Real* Ke, void* ctx);

// ── 弹性力学单元载荷向量 ──
void elasticity_load(const Vec3* coords, std::size_t n_nodes, std::size_t dofs_per_node, Real* Fe, void* ctx);

}  // namespace fem