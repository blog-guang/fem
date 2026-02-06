#include <gtest/gtest.h>
#include "physics/heat_conduction.h"
#include "physics/elasticity.h"
#include "assembly/assembler.h"
#include "assembly/sparse_matrix.h"
#include "assembly/boundary_condition.h"
#include "mesh/mesh_generator.h"
#include "element/element_base.h"
#include "solver/cg.h"
#include "core/types.h"

using namespace fem;

// ── HeatConduction Tests ──
TEST(PhysicsTest, HeatConduction_MatrixSymmetry) {
    // 测试热传导单元刚度矩阵的对称性
    Vec3 coords[3] = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}};
    HeatMaterial mat;
    mat.conductivity = 1.0;
    
    Real Ke[9];
    heat_stiffness(coords, 3, 1, Ke, &mat);
    
    // 检查矩阵是否对称
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            EXPECT_NEAR(Ke[i*3 + j], Ke[j*3 + i], 1e-12);
        }
    }
}

TEST(PhysicsTest, HeatConduction_LoadPositivity) {
    // 测试热传导单元载荷向量的正定性
    Vec3 coords[3] = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}};
    HeatMaterial mat;
    mat.source = 1.0;  // 正的热源
    
    Real Fe[3];
    heat_load(coords, 3, 1, Fe, &mat);
    
    // 检查载荷向量是否为正
    for (int i = 0; i < 3; ++i) {
        EXPECT_GT(Fe[i], 0.0);
    }
}

// ── Elasticity Tests ──
TEST(PhysicsTest, Elasticity_MatrixSymmetry) {
    // 测试弹性力学单元刚度矩阵的对称性
    Vec3 coords[3] = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}};
    ElasticCtx ctx;
    ctx.mat.E = 1e6;
    ctx.mat.poisson = 0.3;
    ctx.mat.thickness = 1.0;
    
    Real Ke[36];  // 6x6 for 3 nodes with 2 DOFs each
    elasticity_stiffness(coords, 3, 2, Ke, &ctx);
    
    // 检查矩阵是否对称
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            EXPECT_NEAR(Ke[i*6 + j], Ke[j*6 + i], 1e-10);
        }
    }
}

TEST(PhysicsTest, Elasticity_RigidBodyModes) {
    // 测试刚体运动模式（理论上应该有3个零特征值）
    // 这个测试较复杂，我们测试刚体位移不会产生内力
    Vec3 coords[3] = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}};
    ElasticCtx ctx;
    ctx.mat.E = 1e6;
    ctx.mat.poisson = 0.3;
    ctx.mat.thickness = 1.0;
    
    Real Ke[36];  // 6x6 for 3 nodes with 2 DOFs each
    elasticity_stiffness(coords, 3, 2, Ke, &ctx);
    
    // 测试刚体位移（平移和旋转）应该不产生力
    // 平移模式
    Real disp_x[6] = {1.0, 0.0, 1.0, 0.0, 1.0, 0.0};  // x方向平移
    Real disp_y[6] = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0};  // y方向平移
    // 刚体旋转（围绕原点，小角度θ=0.01）：u=-y*θ, v=x*θ
    // 节点0(0,0): u=0, v=0; 节点1(1,0): u=0, v=0.01; 节点2(0,1): u=-0.01, v=0
    Real theta = 0.01;
    Real rot[6] = {0.0, 0.0, 0.0, theta, -theta, 0.0};  // 小角度旋转 (绕原点)
    
    // 计算 K*u (应该接近零向量)
    Real force_x[6] = {0.0}, force_y[6] = {0.0}, force_rot[6] = {0.0};
    
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            force_x[i] += Ke[i*6 + j] * disp_x[j];
            force_y[i] += Ke[i*6 + j] * disp_y[j];
            force_rot[i] += Ke[i*6 + j] * rot[j];
        }
    }
    
    // 刚体位移不应产生内力（或非常小的数值误差）
    for (int i = 0; i < 6; ++i) {
        EXPECT_NEAR(force_x[i], 0.0, 1e-10);
        EXPECT_NEAR(force_y[i], 0.0, 1e-10);
        EXPECT_NEAR(force_rot[i], 0.0, 1e-8);  // 旋转模式可能有稍大一些的误差
    }
}

// ── BoundaryCondition Vector Field Tests ──
TEST(PhysicsTest, VectorFieldComponentBC) {
    // 测试矢量场分量边界条件
    Mesh mesh;
    mesh.add_node({0.0, 0.0, 0.0});
    mesh.add_node({1.0, 0.0, 0.0});
    mesh.add_node({0.0, 1.0, 0.0});
    
    // 添加边界标记
    mesh.add_boundary("left", {0});
    mesh.add_boundary("right", {1});
    
    // 创建简单的2x2系统，每节点2个自由度
    CSRMatrix K;
    K.rows = 6;  // 3 nodes * 2 DOFs = 6
    K.row_ptr.resize(7);
    K.col_idx.clear();
    K.values.clear();
    
    // 构建一个简单的对角占优矩阵
    K.row_ptr[0] = 0; K.row_ptr[1] = 1; K.row_ptr[2] = 2; 
    K.row_ptr[3] = 3; K.row_ptr[4] = 4; K.row_ptr[5] = 5; K.row_ptr[6] = 6;
    
    for (int i = 0; i < 6; ++i) {
        K.col_idx.push_back(i);
        K.values.push_back(1.0);  // 对角元素
    }
    
    Vector F(6, 1.0);  // 载荷向量
    
    // 测试固定左侧节点的x分量 (component = 0)
    BoundaryCondition bc_x{BCType::Dirichlet, "left", 5.0, 0};  // 固定x分量为5
    apply_dirichlet(K, F, mesh, bc_x, 2);  // 2 DOFs per node
    
    // 检查第0个自由度 (节点0的x分量) 是否被正确设置
    // 行0应该是 [1, 0, 0, 0, 0, 0]，F[0] 应该是 5.0
    EXPECT_NEAR(F[0], 5.0, 1e-12);
    
    // 重置并测试固定y分量 (component = 1)
    F.assign(6, 1.0);
    BoundaryCondition bc_y{BCType::Dirichlet, "left", 3.0, 1};  // 固定y分量为3
    apply_dirichlet(K, F, mesh, bc_y, 2);  // 2 DOFs per node
    
    // 检查第1个自由度 (节点0的y分量) 是否被正确设置为 3.0
    EXPECT_NEAR(F[1], 3.0, 1e-12);
}

// ── Combined Physics Test: Simple Heat Conduction ──
TEST(PhysicsTest, HeatConduction_SimpleAssembly) {
    // 使用简单网格测试热传导装配
    Mesh mesh = generate_unit_square_tri(2, 2);
    auto elem = create_element(ElementType::Triangle2D);
    Assembler assembler(mesh, *elem, 1);
    
    HeatMaterial mat;
    mat.conductivity = 1.0;
    mat.source = 1.0;
    
    COOMatrix K_coo;
    Vector F;
    assembler.assemble(heat_stiffness, heat_load, K_coo, F, &mat);
    
    // 确保装配成功，矩阵和向量有合理大小
    EXPECT_EQ(F.size(), mesh.num_nodes());
    EXPECT_GT(K_coo.values.size(), 0);
    
    // CSR转换
    CSRMatrix K = coo_to_csr(K_coo);
    EXPECT_EQ(K.rows, mesh.num_nodes());
}

// ── Combined Physics Test: Simple Elasticity ──
TEST(PhysicsTest, Elasticity_SimpleAssembly) {
    // 使用简单网格测试弹性力学装配
    Mesh mesh = generate_unit_square_tri(2, 2);
    auto elem = create_element(ElementType::Triangle2D);
    Assembler assembler(mesh, *elem, 2);  // 2 DOFs per node
    
    ElasticCtx ctx;
    ctx.mat.E = 1e6;
    ctx.mat.poisson = 0.3;
    ctx.mat.thickness = 1.0;
    ctx.load.fx = 0.0;
    ctx.load.fy = -1.0;
    
    COOMatrix K_coo;
    Vector F;
    assembler.assemble(elasticity_stiffness, elasticity_load, K_coo, F, &ctx);
    
    // 确保装配成功，矩阵和向量有合理大小
    EXPECT_EQ(F.size(), mesh.num_nodes() * 2);  // 2 DOFs per node
    EXPECT_GT(K_coo.values.size(), 0);
    
    // CSR转换
    CSRMatrix K = coo_to_csr(K_coo);
    EXPECT_EQ(K.rows, mesh.num_nodes() * 2);  // 2 DOFs per node
}