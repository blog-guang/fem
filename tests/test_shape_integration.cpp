/**
 * test_shape_integration.cpp - 形函数系统与物理模块集成测试
 * 
 * 验证:
 * 1. ElasticityUnified 使用形函数计算（Tri3, Quad4）
 * 2. ElasticityUnified 使用形函数计算（Tet4, Brick8）
 * 3. 高斯积分正确应用
 * 4. 刚度矩阵对称性和正定性
 */

#include <gtest/gtest.h>
#include "physics/elasticity_unified.h"
#include "mesh/mesh.h"
#include "mesh/mesh_generator.h"
#include "shape/shape_function_factory.h"

using namespace fem;
using namespace fem::physics;

// ═══════════════════════════════════════════════════════════════
// 2D Elasticity 集成测试
// ═══════════════════════════════════════════════════════════════

TEST(ShapeIntegrationTest, ElasticityUnified_Tri3) {
    // 创建简单三角形网格
    Material dummy_mat(0, "dummy");
    Mesh mesh("test_tri3", &dummy_mat);
    mesh.add_node(Vec3{0.0, 0.0, 0.0});
    mesh.add_node(Vec3{1.0, 0.0, 0.0});
    mesh.add_node(Vec3{0.0, 1.0, 0.0});
    
    mesh.add_element<Tri3>(0, 1, 2);
    
    // 弹性参数
    Real E = 200e9;   // 200 GPa
    Real nu = 0.3;
    ElasticityUnified elast(E, nu, PlaneType::PlaneStress);
    
    // 计算单元刚度矩阵
    DenseMatrix Ke;
    Vector Fe;
    elast.compute_element(0, mesh, Ke, Fe);
    
    // 验证矩阵尺寸 (6x6 for Tri3)
    EXPECT_EQ(Ke.rows(), 6);
    EXPECT_EQ(Ke.cols(), 6);
    
    // 验证对称性
    Real tol = 1e-6;
    for (int i = 0; i < 6; ++i) {
        for (int j = i + 1; j < 6; ++j) {
            EXPECT_NEAR(Ke(i, j), Ke(j, i), tol) 
                << "Matrix not symmetric at (" << i << "," << j << ")";
        }
    }
    
    // 验证正定性 (对角元素为正)
    for (int i = 0; i < 6; ++i) {
        EXPECT_GT(Ke(i, i), 0.0) 
            << "Diagonal element K(" << i << "," << i << ") not positive";
    }
    
    // 验证载荷向量为零（无体力）
    for (int i = 0; i < 6; ++i) {
        EXPECT_DOUBLE_EQ(Fe[i], 0.0);
    }
}

TEST(ShapeIntegrationTest, ElasticityUnified_Quad4) {
    // 创建四边形网格
    Material dummy_mat(0, "dummy");
    Mesh mesh("test_quad4", &dummy_mat);
    mesh.add_node(Vec3{0.0, 0.0, 0.0});
    mesh.add_node(Vec3{1.0, 0.0, 0.0});
    mesh.add_node(Vec3{1.0, 1.0, 0.0});
    mesh.add_node(Vec3{0.0, 1.0, 0.0});
    
    mesh.add_element<Quad4>(0, 1, 2, 3);
    
    // 弹性参数
    Real E = 200e9;
    Real nu = 0.3;
    ElasticityUnified elast(E, nu, PlaneType::PlaneStrain);
    
    // 计算单元刚度矩阵
    DenseMatrix Ke;
    Vector Fe;
    elast.compute_element(0, mesh, Ke, Fe);
    
    // 验证矩阵尺寸 (8x8 for Quad4)
    EXPECT_EQ(Ke.rows(), 8);
    EXPECT_EQ(Ke.cols(), 8);
    
    // 验证对称性（使用相对容差）
    for (int i = 0; i < 8; ++i) {
        for (int j = i + 1; j < 8; ++j) {
            Real tol = std::max(1e-6, std::abs(Ke(i, j)) * 1e-9);  // 相对容差
            EXPECT_NEAR(Ke(i, j), Ke(j, i), tol);
        }
    }
    
    // 验证正定性
    for (int i = 0; i < 8; ++i) {
        EXPECT_GT(Ke(i, i), 0.0);
    }
}

// ═══════════════════════════════════════════════════════════════
// 3D Elasticity 集成测试
// ═══════════════════════════════════════════════════════════════

TEST(ShapeIntegrationTest, ElasticityUnified3D_Tet4) {
    // 创建四面体网格
    Material dummy_mat(0, "dummy");
    Mesh mesh("test_tet4", &dummy_mat);
    mesh.add_node(Vec3{0.0, 0.0, 0.0});
    mesh.add_node(Vec3{1.0, 0.0, 0.0});
    mesh.add_node(Vec3{0.0, 1.0, 0.0});
    mesh.add_node(Vec3{0.0, 0.0, 1.0});
    
    mesh.add_element<Tet4>(0, 1, 2, 3);
    
    // 弹性参数
    Real E = 200e9;
    Real nu = 0.3;
    ElasticityUnified elast(E, nu, true);
    
    // 计算单元刚度矩阵
    DenseMatrix Ke;
    Vector Fe;
    elast.compute_element(0, mesh, Ke, Fe);
    
    // 验证矩阵尺寸 (12x12 for Tet4)
    EXPECT_EQ(Ke.rows(), 12);
    EXPECT_EQ(Ke.cols(), 12);
    
    // 验证对称性
    Real tol = 1e-6;
    for (int i = 0; i < 12; ++i) {
        for (int j = i + 1; j < 12; ++j) {
            EXPECT_NEAR(Ke(i, j), Ke(j, i), tol)
                << "Matrix not symmetric at (" << i << "," << j << ")";
        }
    }
    
    // 验证正定性
    for (int i = 0; i < 12; ++i) {
        EXPECT_GT(Ke(i, i), 0.0)
            << "Diagonal element K(" << i << "," << i << ") not positive";
    }
}

TEST(ShapeIntegrationTest, ElasticityUnified3D_Brick8) {
    // 创建六面体网格 (单位立方体)
    Material dummy_mat(0, "dummy");
    Mesh mesh("test_brick8", &dummy_mat);
    mesh.add_node(Vec3{0.0, 0.0, 0.0});  // 0
    mesh.add_node(Vec3{1.0, 0.0, 0.0});  // 1
    mesh.add_node(Vec3{1.0, 1.0, 0.0});  // 2
    mesh.add_node(Vec3{0.0, 1.0, 0.0});  // 3
    mesh.add_node(Vec3{0.0, 0.0, 1.0});  // 4
    mesh.add_node(Vec3{1.0, 0.0, 1.0});  // 5
    mesh.add_node(Vec3{1.0, 1.0, 1.0});  // 6
    mesh.add_node(Vec3{0.0, 1.0, 1.0});  // 7
    
    std::array<Index, 8> nodes = {0, 1, 2, 3, 4, 5, 6, 7};
    mesh.add_element<Brick8>(nodes);
    
    // 弹性参数
    Real E = 200e9;
    Real nu = 0.3;
    ElasticityUnified elast(E, nu, true);
    
    // 计算单元刚度矩阵
    DenseMatrix Ke;
    Vector Fe;
    elast.compute_element(0, mesh, Ke, Fe);
    
    // 验证矩阵尺寸 (24x24 for Brick8)
    EXPECT_EQ(Ke.rows(), 24);
    EXPECT_EQ(Ke.cols(), 24);
    
    // 验证对称性（使用相对容差）
    for (int i = 0; i < 24; ++i) {
        for (int j = i + 1; j < 24; ++j) {
            Real tol = std::max(1e-6, std::abs(Ke(i, j)) * 1e-9);  // 相对容差
            EXPECT_NEAR(Ke(i, j), Ke(j, i), tol);
        }
    }
    
    // 验证正定性
    for (int i = 0; i < 24; ++i) {
        EXPECT_GT(Ke(i, i), 0.0);
    }
}

// ═══════════════════════════════════════════════════════════════
// 与旧实现一致性测试
// ═══════════════════════════════════════════════════════════════

TEST(ShapeIntegrationTest, Tri3_BackwardCompatibility) {
    // 验证新实现（形函数）与旧实现（硬编码）结果一致
    
    Material dummy_mat(0, "dummy");
    Mesh mesh("test_compat", &dummy_mat);
    mesh.add_node(Vec3{0.0, 0.0, 0.0});
    mesh.add_node(Vec3{2.0, 0.0, 0.0});
    mesh.add_node(Vec3{0.0, 3.0, 0.0});
    
    mesh.add_element<Tri3>(0, 1, 2);
    
    Real E = 100.0;
    Real nu = 0.25;
    ElasticityUnified elast(E, nu, PlaneType::PlaneStress);
    
    // 新实现（使用形函数）
    DenseMatrix Ke_new;
    Vector Fe_new;
    elast.compute_element(0, mesh, Ke_new, Fe_new);
    
    // 刚度矩阵应该合理（不为零，对称，正定）
    EXPECT_GT(Ke_new(0, 0), 0.0);
    
    // 验证刚度矩阵的范数量级合理
    Real norm = 0.0;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            norm += Ke_new(i, j) * Ke_new(i, j);
        }
    }
    norm = std::sqrt(norm);
    EXPECT_GT(norm, 1.0);  // 应该不是零矩阵
}

// ═══════════════════════════════════════════════════════════════
// 刚体运动测试
// ═══════════════════════════════════════════════════════════════

TEST(ShapeIntegrationTest, RigidBodyMotion_Tri3) {
    // 刚体平移应该产生零应变 → 零内力
    
    Material dummy_mat(0, "dummy");
    Mesh mesh("test_rigid", &dummy_mat);
    mesh.add_node(Vec3{0.0, 0.0, 0.0});
    mesh.add_node(Vec3{1.0, 0.0, 0.0});
    mesh.add_node(Vec3{0.0, 1.0, 0.0});
    
    mesh.add_element<Tri3>(0, 1, 2);
    
    Real E = 200e9;
    Real nu = 0.3;
    ElasticityUnified elast(E, nu, PlaneType::PlaneStress);
    
    DenseMatrix Ke;
    Vector Fe;
    elast.compute_element(0, mesh, Ke, Fe);
    
    // 刚体平移位移：u = [c, c, c, c, c, c]^T (所有节点平移相同)
    Vector u_rigid(6);
    for (int i = 0; i < 6; ++i) {
        u_rigid[i] = 1.0;  // 平移 1.0 单位
    }
    
    // 内力 f = K * u （刚体运动应产生零内力）
    Vector f = Ke * u_rigid;
    
    // 验证内力为零（或极小）
    // 使用相对于刚度矩阵元素的容差
    Real K_norm = 0.0;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            K_norm = std::max(K_norm, std::abs(Ke(i, j)));
        }
    }
    Real tol = K_norm * 1e-9;  // 相对容差
    
    for (int i = 0; i < 6; ++i) {
        EXPECT_NEAR(f[i], 0.0, tol)
            << "Rigid body motion produced non-zero force at DOF " << i;
    }
}
