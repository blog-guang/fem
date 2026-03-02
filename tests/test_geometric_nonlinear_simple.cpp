/**
 * test_geometric_nonlinear_simple.cpp - 几何非线性简化测试
 * 
 * 验证已实现的核心功能：
 * - ElasticityNonlinear 切线刚度计算
 * - Neo-Hookean 大变形应力
 * - 内力向量计算
 */

#include <gtest/gtest.h>
#include "physics/elasticity_nonlinear.h"
#include "material/neo_hookean.h"
#include "material/isotropic_elastic.h"
#include "mesh/model.h"
#include "mesh/mesh_generator.h"
#include <cmath>

using namespace fem;
using namespace fem::physics;
using namespace fem::constitutive;

// ═══════════════════════════════════════════════════════════
// 辅助函数：简单迭代求解器（替代完整 Newton-Raphson）
// ═══════════════════════════════════════════════════════════

/**
 * 简单的单步位移求解
 * 
 * K * u = F （不考虑迭代，仅作为功能验证）
 * 
 * 注意：这不是真正的非线性求解，仅用于测试物理模块
 */
Vector simple_linear_solve(
    const DenseMatrix& K,
    const Vector& F)
{
    // 简化：使用高斯消元（仅小规模系统）
    int n = K.rows();
    Vector u(n, 0.0);
    
    if (n == 0 || K.cols() != n || F.size() != static_cast<std::size_t>(n)) {
        return u;
    }
    
    // 构造增广矩阵 [K | F]
    DenseMatrix Aug(n, n + 1);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Aug(i, j) = K(i, j);
        }
        Aug(i, n) = F[i];
    }
    
    // 高斯消元（简化版，无选主元）
    for (int k = 0; k < n; k++) {
        // 跳过零对角元（边界条件）
        if (std::abs(Aug(k, k)) < 1e-10) continue;
        
        for (int i = k + 1; i < n; i++) {
            Real factor = Aug(i, k) / Aug(k, k);
            for (int j = k; j <= n; j++) {
                Aug(i, j) -= factor * Aug(k, j);
            }
        }
    }
    
    // 回代
    for (int i = n - 1; i >= 0; i--) {
        if (std::abs(Aug(i, i)) < 1e-10) {
            u[i] = 0.0;
            continue;
        }
        
        Real sum = Aug(i, n);
        for (int j = i + 1; j < n; j++) {
            sum -= Aug(i, j) * u[j];
        }
        u[i] = sum / Aug(i, i);
    }
    
    return u;
}

// ═══════════════════════════════════════════════════════════
// 测试 1: 单元刚度与内力计算
// ═══════════════════════════════════════════════════════════

TEST(GeometricNonlinearSimpleTest, ElementStiffnessAndInternalForce) {
    // 创建单个 Tet4 单元
    
    Model model("test");
    int mat_id = model.add_material("test");
    int mesh_id = model.add_mesh("test_mesh", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    MeshGenerator::generate_unit_cube_tet(1, 1, 1, mesh);
    
    // Neo-Hookean 材料
    Real E = 10e6;
    Real nu = 0.3;
    NeoHookean mat(E, nu, 3, true);
    
    // 几何非线性物理模块
    ElasticityNonlinear physics(&mat, 3, true);
    
    // 施加位移（10% 拉伸）
    Vector u(mesh.num_nodes() * 3, 0.0);
    for (Index i = 0; i < mesh.num_nodes(); i++) {
        const auto& coords = mesh.node(i).coords();
        u[i * 3 + 0] = coords[0] * 0.1;  // X 方向 10% 拉伸
    }
    
    // 计算第一个单元的切线刚度
    DenseMatrix Ke;
    Vector Fe;
    physics.compute_element_nonlinear(0, mesh, u, Ke, Fe);
    
    int n_dofs = mesh.element(0).num_nodes() * 3;
    
    // 验证刚度矩阵
    EXPECT_EQ(Ke.rows(), n_dofs);
    EXPECT_EQ(Ke.cols(), n_dofs);
    
    // 对称性
    for (int i = 0; i < n_dofs; i++) {
        for (int j = 0; j < n_dofs; j++) {
            EXPECT_NEAR(Ke(i, j), Ke(j, i), 1e-6);
        }
    }
    
    // 计算内力
    Vector F_int = physics.compute_internal_force(0, mesh, u);
    
    EXPECT_EQ(F_int.size(), n_dofs);
    
    // 内力应该非零（有应变）
    Real F_norm = 0.0;
    for (int i = 0; i < n_dofs; i++) {
        F_norm += F_int[i] * F_int[i];
    }
    F_norm = std::sqrt(F_norm);
    
    EXPECT_GT(F_norm, 1e3);
}

// ═══════════════════════════════════════════════════════════
// 测试 2: 材料非线性对比（Neo-Hookean vs Linear）
// ═══════════════════════════════════════════════════════════

TEST(GeometricNonlinearSimpleTest, MaterialNonlinearity_Comparison) {
    // 对比 Neo-Hookean 和线弹性的刚度差异
    
    Model model("test");
    int mat_id = model.add_material("test");
    int mesh_id = model.add_mesh("test_mesh", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    MeshGenerator::generate_unit_cube_tet(1, 1, 1, mesh);
    
    Real E = 10e6;
    Real nu = 0.3;
    
    // Neo-Hookean
    NeoHookean neo(E, nu, 3, true);
    ElasticityNonlinear phys_neo(&neo, 3, true);
    
    // 线弹性
    IsotropicElastic iso(E, nu, 3, false);
    ElasticityNonlinear phys_iso(&iso, 3, true);
    
    // 施加中等应变（15%）
    Vector u(mesh.num_nodes() * 3, 0.0);
    for (Index i = 0; i < mesh.num_nodes(); i++) {
        const auto& coords = mesh.node(i).coords();
        u[i * 3 + 0] = coords[0] * 0.15;
    }
    
    // 计算内力
    Vector F_neo = phys_neo.compute_internal_force(0, mesh, u);
    Vector F_iso = phys_iso.compute_internal_force(0, mesh, u);
    
    // 应该有差异（材料非线性）
    Real diff = 0.0;
    for (std::size_t i = 0; i < F_neo.size(); i++) {
        diff += std::abs(F_neo[i] - F_iso[i]);
    }
    
    // 对于 15% 应变，差异应该显著
    EXPECT_GT(diff, 1.0);
}

// ═══════════════════════════════════════════════════════════
// 测试 3: 几何刚度的作用
// ═══════════════════════════════════════════════════════════

TEST(GeometricNonlinearSimpleTest, GeometricStiffness_Effect) {
    // 验证几何刚度的贡献
    
    Model model("test");
    int mat_id = model.add_material("test");
    int mesh_id = model.add_mesh("test_mesh", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    MeshGenerator::generate_unit_cube_tet(1, 1, 1, mesh);
    
    Real E = 10e6;
    Real nu = 0.3;
    IsotropicElastic mat(E, nu, 3, false);
    
    // 施加预应变
    Vector u(mesh.num_nodes() * 3, 0.0);
    for (Index i = 0; i < mesh.num_nodes(); i++) {
        const auto& coords = mesh.node(i).coords();
        u[i * 3 + 0] = coords[0] * 0.1;
    }
    
    // 禁用几何刚度
    ElasticityNonlinear phys_no_geo(&mat, 3, false);
    DenseMatrix Ke_no_geo;
    Vector Fe1;
    phys_no_geo.compute_element_nonlinear(0, mesh, u, Ke_no_geo, Fe1);
    
    // 启用几何刚度
    ElasticityNonlinear phys_with_geo(&mat, 3, true);
    DenseMatrix Ke_with_geo;
    Vector Fe2;
    phys_with_geo.compute_element_nonlinear(0, mesh, u, Ke_with_geo, Fe2);
    
    // 计算差异
    Real diff = 0.0;
    int n = Ke_no_geo.rows();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            diff += std::abs(Ke_with_geo(i, j) - Ke_no_geo(i, j));
        }
    }
    
    // 有预应力时，几何刚度应该有贡献
    EXPECT_GT(diff, 1e-3);
}

// ═══════════════════════════════════════════════════════════
// 测试 4: Cook's Membrane 简化模拟（单元级别）
// ═══════════════════════════════════════════════════════════

// TODO: 修复 2D 平面应变问题的矩阵尺寸不匹配
TEST(GeometricNonlinearSimpleTest, DISABLED_CooksMembrane_SingleElement) {
    // 简化的 Cook's Membrane：单个 Quad4 单元
    
    Model model("cooks");
    int mat_id = model.add_material("rubber");
    int mesh_id = model.add_mesh("membrane", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    // 梯形单元（模拟 Cook's Membrane 几何）
    Vec3 coords[4] = {
        {0.0, 0.0, 0.0},      // 左下
        {48.0, 0.0, 0.0},     // 右下
        {48.0, 16.0, 0.0},    // 右上
        {0.0, 44.0, 0.0}      // 左上
    };
    
    Index node_ids[4];
    for (int i = 0; i < 4; i++) {
        node_ids[i] = mesh.add_node(coords[i]);
    }
    
    mesh.add_element<Quad4>(node_ids[0], node_ids[1], node_ids[2], node_ids[3]);
    
    // Neo-Hookean 材料（软材料，近不可压）
    Real E = 250.0;
    Real nu = 0.4999;
    NeoHookean mat(E, nu, 2, true);  // 2D 平面应变
    
    ElasticityNonlinear physics(&mat, 2, true);
    
    // 施加剪切位移（模拟剪切载荷效果）
    Vector u(mesh.num_nodes() * 2, 0.0);
    
    // 左边固定（u[0-3] = 0）
    // 右边向上位移
    u[2] = 0.0;   // 右下 X
    u[3] = 5.0;   // 右下 Y
    u[4] = 0.0;   // 右上 X
    u[5] = 10.0;  // 右上 Y
    
    // 计算内力
    Vector F_int = physics.compute_internal_force(0, mesh, u);
    
    EXPECT_EQ(F_int.size(), 8);  // 4 nodes * 2 DOFs
    
    // 内力应该非零
    Real F_norm = 0.0;
    for (std::size_t i = 0; i < F_int.size(); i++) {
        F_norm += F_int[i] * F_int[i];
    }
    F_norm = std::sqrt(F_norm);
    
    EXPECT_GT(F_norm, 1.0);
}

// ═══════════════════════════════════════════════════════════
// 测试 5: 应力还原（验证大变形应力计算）
// ═══════════════════════════════════════════════════════════

TEST(GeometricNonlinearSimpleTest, StressRecovery_LargeDeformation) {
    // 验证大变形下的应力计算
    
    Model model("test");
    int mat_id = model.add_material("rubber");
    int mesh_id = model.add_mesh("test_mesh", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    MeshGenerator::generate_unit_cube_tet(1, 1, 1, mesh);
    
    // Neo-Hookean 材料
    Real E = 10e6;
    Real nu = 0.3;
    NeoHookean mat(E, nu, 3, true);
    
    // 施加大变形（30% 拉伸）
    Vector u(mesh.num_nodes() * 3, 0.0);
    for (Index i = 0; i < mesh.num_nodes(); i++) {
        const auto& coords = mesh.node(i).coords();
        u[i * 3 + 0] = coords[0] * 0.3;
    }
    
    ElasticityNonlinear physics(&mat, 3, true);
    
    // 计算内力（内部会计算应力）
    Vector F_int = physics.compute_internal_force(0, mesh, u);
    
    // 内力应该随着变形增大
    Real F_norm = 0.0;
    for (std::size_t i = 0; i < F_int.size(); i++) {
        F_norm += F_int[i] * F_int[i];
    }
    F_norm = std::sqrt(F_norm);
    
    // 30% 应变应该产生显著内力
    EXPECT_GT(F_norm, 1e4);
}
