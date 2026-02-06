#include <gtest/gtest.h>
#include "assembly/assembler.h"
#include "mesh/model.h"
#include "mesh/mesh_generator.h"

using namespace fem;

// ═══ Assembler Tests ═══

TEST(AssemblerTest, Construction) {
    Model model("test");
    int mat_id = model.add_material("mat");
    int mesh_id = model.add_mesh("mesh", mat_id);
    
    Mesh& mesh = model.mesh(mesh_id);
    MeshGenerator::generate_unit_square_tri(2, 2, mesh);
    
    // 标量场
    Assembler assembler1(model, 1);
    EXPECT_EQ(assembler1.dofs_per_node(), 1u);
    EXPECT_EQ(assembler1.num_dofs(), mesh.num_nodes());
    
    // 2D 矢量场
    Assembler assembler2(model, 2);
    EXPECT_EQ(assembler2.dofs_per_node(), 2u);
    EXPECT_EQ(assembler2.num_dofs(), mesh.num_nodes() * 2);
}

TEST(AssemblerTest, AssembleScalar) {
    // 创建简单的 2x2 网格
    Model model("test");
    int mat_id = model.add_material("mat");
    int mesh_id = model.add_mesh("mesh", mat_id);
    
    Mesh& mesh = model.mesh(mesh_id);
    MeshGenerator::generate_unit_square_tri(2, 2, mesh);
    
    Assembler assembler(model, 1);
    
    // 简单的单元矩阵: Ke = [[1,0,0],[0,1,0],[0,0,1]] (单位矩阵)
    // 简单的单元向量: Fe = [1,1,1]
    auto elem_func = [](Index elem_id, const Mesh& m, DenseMatrix& Ke, Vector& Fe) {
        Ke.identity();  // 单位矩阵
        for (std::size_t i = 0; i < Fe.size(); ++i) {
            Fe[i] = 1.0;  // 单位载荷
        }
    };
    
    assembler.assemble(elem_func);
    
    auto K = assembler.matrix();
    const auto& F = assembler.rhs();
    
    EXPECT_GT(K.nnz(), 0u);
    EXPECT_EQ(F.size(), mesh.num_nodes());
    
    // 检查载荷向量非零
    Real sum_F = 0.0;
    for (std::size_t i = 0; i < F.size(); ++i) {
        sum_F += F[i];
    }
    EXPECT_GT(sum_F, 0.0);
}

TEST(AssemblerTest, AssembleVector) {
    // 创建 2x2 网格
    Model model("test");
    int mat_id = model.add_material("mat");
    int mesh_id = model.add_mesh("mesh", mat_id);
    
    Mesh& mesh = model.mesh(mesh_id);
    MeshGenerator::generate_unit_square_tri(2, 2, mesh);
    
    // 2D 矢量场 (每个节点2个自由度)
    Assembler assembler(model, 2);
    
    auto elem_func = [](Index elem_id, const Mesh& m, DenseMatrix& Ke, Vector& Fe) {
        // 6x6 单位矩阵 (3个节点 * 2个自由度)
        Ke.identity();
        for (std::size_t i = 0; i < Fe.size(); ++i) {
            Fe[i] = 1.0;
        }
    };
    
    assembler.assemble(elem_func);
    
    auto K = assembler.matrix();
    const auto& F = assembler.rhs();
    
    EXPECT_EQ(assembler.num_dofs(), mesh.num_nodes() * 2);
    EXPECT_EQ(F.size(), mesh.num_nodes() * 2);
    EXPECT_GT(K.nnz(), 0u);
}

TEST(AssemblerTest, ApplyDirichlet) {
    // 创建 2x2 网格
    Model model("test");
    int mat_id = model.add_material("mat");
    int mesh_id = model.add_mesh("mesh", mat_id);
    
    Mesh& mesh = model.mesh(mesh_id);
    MeshGenerator::generate_unit_square_tri(2, 2, mesh);
    MeshGenerator::identify_boundaries_2d(mesh);
    
    Assembler assembler(model, 1);
    
    // 装配简单系统
    auto elem_func = [](Index elem_id, const Mesh& m, DenseMatrix& Ke, Vector& Fe) {
        Ke.identity();
        for (std::size_t i = 0; i < Fe.size(); ++i) {
            Fe[i] = 1.0;
        }
    };
    
    assembler.assemble(elem_func);
    
    // 应用 Dirichlet 边界条件: left boundary, dof=0, value=5.0
    std::vector<DirichletBC> bcs;
    bcs.push_back({"left", 0, 5.0});
    
    assembler.apply_dirichlet(bcs);
    
    const auto& F = assembler.rhs();
    const auto& left_nodes = mesh.boundary("left");
    
    // 检查左边界节点的载荷值被设置为5.0
    for (Index node_id : left_nodes) {
        EXPECT_DOUBLE_EQ(F[node_id], 5.0);
    }
}

TEST(AssemblerTest, ApplyDirichletVector) {
    // 创建 2x2 网格
    Model model("test");
    int mat_id = model.add_material("mat");
    int mesh_id = model.add_mesh("mesh", mat_id);
    
    Mesh& mesh = model.mesh(mesh_id);
    MeshGenerator::generate_unit_square_tri(2, 2, mesh);
    MeshGenerator::identify_boundaries_2d(mesh);
    
    // 2D 矢量场
    Assembler assembler(model, 2);
    
    auto elem_func = [](Index elem_id, const Mesh& m, DenseMatrix& Ke, Vector& Fe) {
        Ke.identity();
        for (std::size_t i = 0; i < Fe.size(); ++i) {
            Fe[i] = 1.0;
        }
    };
    
    assembler.assemble(elem_func);
    
    // 应用 Dirichlet 边界条件:
    // - left boundary, dof=0 (x方向), value=1.0
    // - right boundary, dof=1 (y方向), value=2.0
    std::vector<DirichletBC> bcs;
    bcs.push_back({"left", 0, 1.0});
    bcs.push_back({"right", 1, 2.0});
    
    assembler.apply_dirichlet(bcs);
    
    const auto& F = assembler.rhs();
    const auto& left_nodes = mesh.boundary("left");
    const auto& right_nodes = mesh.boundary("right");
    
    // 检查左边界节点的 x 自由度 (dof=0) 被设置为 1.0
    for (Index node_id : left_nodes) {
        Index dof_x = node_id * 2 + 0;
        EXPECT_DOUBLE_EQ(F[dof_x], 1.0);
    }
    
    // 检查右边界节点的 y 自由度 (dof=1) 被设置为 2.0
    for (Index node_id : right_nodes) {
        Index dof_y = node_id * 2 + 1;
        EXPECT_DOUBLE_EQ(F[dof_y], 2.0);
    }
}

TEST(AssemblerTest, Clear) {
    Model model("test");
    int mat_id = model.add_material("mat");
    int mesh_id = model.add_mesh("mesh", mat_id);
    
    Mesh& mesh = model.mesh(mesh_id);
    MeshGenerator::generate_unit_square_tri(2, 2, mesh);
    
    Assembler assembler(model, 1);
    
    auto elem_func = [](Index elem_id, const Mesh& m, DenseMatrix& Ke, Vector& Fe) {
        Ke.identity();
        for (std::size_t i = 0; i < Fe.size(); ++i) {
            Fe[i] = 1.0;
        }
    };
    
    assembler.assemble(elem_func);
    EXPECT_GT(assembler.matrix().nnz(), 0u);
    
    assembler.clear();
    
    // 清空后再次装配
    assembler.assemble(elem_func);
    EXPECT_GT(assembler.matrix().nnz(), 0u);
}
