/**
 * dof_handler_demo.cpp - 演示新的 DOF 管理架构
 * 
 * 展示:
 * 1. DofHandler DOF 编号
 * 2. SparseMatrixPattern 预分配
 * 3. 基于 Pattern 的矩阵初始化
 */

#include "mesh/model.h"
#include "mesh/mesh_generator.h"
#include "mesh/dof_handler.h"
#include "math/sparse_matrix.h"
#include "core/logger.h"
#include <iostream>

using namespace fem;

int main() {
    FEM_INFO("=== DofHandler Architecture Demo ===\n");
    
    // 1. 创建模型
    Model model("demo");
    model.add_material("steel");
    model.add_mesh("beam", 0);
    
    // 2. 生成简单网格 (3x3 四边形网格)
    auto& mesh = model.mesh(0);
    MeshGenerator::generate_unit_square_quad(3, 3, mesh);
    
    FEM_INFO("Mesh: " + std::to_string(mesh.num_nodes()) + " nodes, " +
             std::to_string(mesh.num_elements()) + " elements");
    
    // 3. 创建 DOF 处理器
    DofHandler dof_handler(model);
    
    // 3.1 设置有限元空间 (3D 位移场)
    FiniteElementSpace fe_space(FieldType::Vector3D, 1, 3);
    dof_handler.distribute_dofs(fe_space);
    
    FEM_INFO("DOFs: " + std::to_string(dof_handler.n_dofs()) + 
             ", DOFs/node: " + std::to_string(dof_handler.dofs_per_node()));
    
    // 4. 构建稀疏模式
    FEM_INFO("\nBuilding sparsity pattern...");
    auto pattern = dof_handler.make_sparsity_pattern();
    
    pattern.print_info();
    
    // 5. 基于 Pattern 创建矩阵（预分配内存）
    FEM_INFO("\nCreating matrix from pattern (pre-allocated)...");
    auto K = pattern.create_matrix();
    
    K.print("Stiffness Matrix");
    
    // 6. 添加 Dirichlet 边界条件
    FEM_INFO("\nAdding Dirichlet BCs...");
    // 固定左边界 (x≈0) 的所有节点
    for (Index node_id = 0; node_id < mesh.num_nodes(); node_id++) {
        const auto& node = mesh.node(node_id);
        auto coords = node.coords();
        if (std::abs(coords[0]) < 1e-6) {  // x 坐标
            // 固定 u, v, w
            dof_handler.add_dirichlet_bc(node_id, 0, 0.0);  // u = 0
            dof_handler.add_dirichlet_bc(node_id, 1, 0.0);  // v = 0
            dof_handler.add_dirichlet_bc(node_id, 2, 0.0);  // w = 0
        }
    }
    
    FEM_INFO("Dirichlet DOFs: " + std::to_string(dof_handler.dirichlet_dofs().size()));
    
    // 7. 演示节点 DOF 映射
    FEM_INFO("\nNode DOF mapping example:");
    if (mesh.num_nodes() > 0) {
        auto dofs = dof_handler.node_dofs(0);
        FEM_INFO("  Node 0 DOFs: ");
        for (size_t i = 0; i < dofs.size(); i++) {
            std::cout << "    component " << i << " -> global DOF " << dofs[i] << std::endl;
        }
    }
    
    // 8. 演示单元 DOF 映射
    FEM_INFO("\nElement DOF mapping example:");
    if (mesh.num_elements() > 0) {
        auto dofs = dof_handler.element_dofs(0);
        std::cout << "  Element 0 has " << dofs.size() << " DOFs:" << std::endl;
        for (size_t i = 0; i < std::min(dofs.size(), size_t(5)); i++) {
            std::cout << "    local " << i << " -> global DOF " << dofs[i] << std::endl;
        }
        if (dofs.size() > 5) {
            std::cout << "    ... and " << (dofs.size() - 5) << " more" << std::endl;
        }
    }
    
    FEM_INFO("\n=== Demo Complete ===");
    
    return 0;
}
