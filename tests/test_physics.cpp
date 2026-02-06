#include <gtest/gtest.h>
#include "physics/heat.h"
#include "physics/elasticity_v2.h"
#include "mesh/model.h"
#include "mesh/mesh_generator.h"

using namespace fem;
using namespace fem::physics;

// ═══ HeatConduction Tests ═══
TEST(HeatConductionTest, Construction) {
    HeatConduction heat(1.5, 10.0);
    
    EXPECT_DOUBLE_EQ(heat.conductivity(), 1.5);
    EXPECT_DOUBLE_EQ(heat.source(), 10.0);
}

TEST(HeatConductionTest, SimpleElement) {
    Model model("test");
    int mat_id = model.add_material("mat");
    int mesh_id = model.add_mesh("mesh", mat_id);
    
    Mesh& mesh = model.mesh(mesh_id);
    MeshGenerator::generate_unit_square_tri(2, 2, mesh);
    
    HeatConduction heat(1.0, 1.0);
    
    DenseMatrix Ke(3, 3);
    Vector Fe(3);
    
    // 计算第一个单元
    heat.compute_element(0, mesh, Ke, Fe);
    
    // 刚度矩阵应该是对称的
    EXPECT_NEAR(Ke(0, 1), Ke(1, 0), 1e-12);
    EXPECT_NEAR(Ke(0, 2), Ke(2, 0), 1e-12);
    EXPECT_NEAR(Ke(1, 2), Ke(2, 1), 1e-12);
    
    // 对角元应该是正的
    EXPECT_GT(Ke(0, 0), 0.0);
    EXPECT_GT(Ke(1, 1), 0.0);
    EXPECT_GT(Ke(2, 2), 0.0);
    
    // 载荷向量应该是正的（源项为1）
    EXPECT_GT(Fe[0], 0.0);
    EXPECT_GT(Fe[1], 0.0);
    EXPECT_GT(Fe[2], 0.0);
}

// ═══ Elasticity2D Tests ═══
TEST(Elasticity2DTest, Construction) {
    Elasticity2D elast(1000.0, 0.3, PlaneType::PlaneStress);
    
    EXPECT_DOUBLE_EQ(elast.youngs_modulus(), 1000.0);
    EXPECT_DOUBLE_EQ(elast.poissons_ratio(), 0.3);
    EXPECT_EQ(elast.plane_type(), PlaneType::PlaneStress);
}

TEST(Elasticity2DTest, SimpleElement) {
    Model model("test");
    int mat_id = model.add_material("mat");
    int mesh_id = model.add_mesh("mesh", mat_id);
    
    Mesh& mesh = model.mesh(mesh_id);
    MeshGenerator::generate_unit_square_tri(2, 2, mesh);
    
    Elasticity2D elast(1000.0, 0.3, PlaneType::PlaneStress);
    
    DenseMatrix Ke(6, 6);  // 3 节点 * 2 DOF
    Vector Fe(6);
    
    // 计算第一个单元
    elast.compute_element(0, mesh, Ke, Fe);
    
    // 刚度矩阵应该是对称的
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            EXPECT_NEAR(Ke(i, j), Ke(j, i), 1e-10) << "i=" << i << ", j=" << j;
        }
    }
    
    // 对角元应该是正的
    for (int i = 0; i < 6; ++i) {
        EXPECT_GT(Ke(i, i), 0.0) << "Diagonal element " << i << " should be positive";
    }
    
    // 载荷向量应该是0（无体力）
    for (int i = 0; i < 6; ++i) {
        EXPECT_DOUBLE_EQ(Fe[i], 0.0);
    }
}

TEST(Elasticity2DTest, MatrixPositiveDefinite) {
    Model model("test");
    int mat_id = model.add_material("mat");
    int mesh_id = model.add_mesh("mesh", mat_id);
    
    Mesh& mesh = model.mesh(mesh_id);
    MeshGenerator::generate_unit_square_tri(2, 2, mesh);
    
    Elasticity2D elast(1000.0, 0.3, PlaneType::PlaneStress);
    
    DenseMatrix Ke(6, 6);
    Vector Fe(6);
    
    elast.compute_element(0, mesh, Ke, Fe);
    
    // 检查矩阵的正定性：所有主对角元素应该是正的
    for (int i = 0; i < 6; ++i) {
        EXPECT_GT(Ke(i, i), 0.0);
    }
    
    // 简单的正定性检查：对于任意非零向量 v，v^T K v > 0
    Vector v(6);
    for (int i = 0; i < 6; ++i) {
        v[i] = 1.0;
    }
    
    Real result = 0.0;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            result += v[i] * Ke(i, j) * v[j];
        }
    }
    
    EXPECT_GT(result, 0.0) << "Matrix should be positive definite";
}
