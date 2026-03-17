/**
 * dof_handler_benchmark.cpp - DofHandler 性能测试
 * 
 * 测试:
 * 1. 初始化时间
 * 2. 稀疏模式构建时间
 * 3. DOF 查询性能
 */

#include "mesh/model.h"
#include "mesh/mesh_generator.h"
#include "mesh/dof_handler.h"
#include "core/logger.h"
#include "core/timer.h"
#include <iostream>

using namespace fem;

struct TestConfig {
    std::string name;
    int nx, ny, nz;
};

void run_benchmark(const TestConfig& config) {
    FEM_INFO("\n" + std::string(60, '='));
    FEM_INFO("测试网格: " + config.name);
    FEM_INFO(std::string(60, '='));
    
    // 创建模型
    Model model("benchmark");
    model.add_material("steel");
    model.add_mesh("domain", 0);
    
    auto& mesh = model.mesh(0);
    
    if (config.nz > 0) {
        MeshGenerator::generate_unit_cube_tet(config.nx, config.ny, config.nz, mesh);
    } else {
        MeshGenerator::generate_unit_square_quad(config.nx, config.ny, mesh);
    }
    
    FEM_INFO("网格: " + std::to_string(mesh.num_nodes()) + " 节点, " +
             std::to_string(mesh.num_elements()) + " 单元");
    
    Timer timer;
    
    // 1. 创建 DofHandler
    timer.start();
    DofHandler dof_handler(model);
    FiniteElementSpace fe_space(FieldType::Vector3D, 1);
    dof_handler.distribute_dofs(fe_space);
    auto t_init = timer.elapsed_ms();
    
    FEM_INFO("初始化: " + std::to_string(t_init) + " ms");
    FEM_INFO("  总 DOFs: " + std::to_string(dof_handler.n_dofs()));
    
    // 2. 构建稀疏模式
    timer.start();
    auto pattern = dof_handler.make_sparsity_pattern();
    auto t_pattern = timer.elapsed_ms();
    
    FEM_INFO("稀疏模式构建: " + std::to_string(t_pattern) + " ms");
    pattern.print_info();
    
    // 3. 测试 DOF 查询性能
    timer.start();
    const size_t num_queries = 10000;
    size_t sum = 0;
    
    for (size_t i = 0; i < num_queries; i++) {
        Index node_id = i % mesh.num_nodes();
        for (int comp = 0; comp < dof_handler.dofs_per_node(); comp++) {
            sum += dof_handler.node_dof(node_id, comp);
        }
    }
    
    auto t_query = timer.elapsed_ms();
    
    FEM_INFO("DOF 查询 (" + std::to_string(num_queries) + " 次): " + 
             std::to_string(t_query) + " ms");
    FEM_INFO("  查询速率: " + 
             std::to_string(num_queries / (t_query / 1000.0)) + " queries/s");
    
    // 4. 测试单元 DOF 查询
    timer.start();
    std::vector<Index> elem_dofs;
    for (Index elem_id = 0; elem_id < mesh.num_elements(); elem_id++) {
        dof_handler.element_dofs(elem_id, elem_dofs);
        sum += elem_dofs.size();
    }
    auto t_elem = timer.elapsed_ms();
    
    FEM_INFO("单元 DOF 查询 (" + std::to_string(mesh.num_elements()) + " 单元): " + 
             std::to_string(t_elem) + " ms");
}

int main() {
    FEM_INFO("=== DofHandler 性能测试 ===\n");
    
    std::vector<TestConfig> configs = {
        {"小网格 (10x10)", 10, 10, 0},
        {"中等网格 (30x30)", 30, 30, 0},
        {"大网格 (50x50)", 50, 50, 0},
        {"3D 网格 (10x10x10)", 10, 10, 10}
    };
    
    for (const auto& config : configs) {
        run_benchmark(config);
    }
    
    FEM_INFO("\n=== Benchmark 完成 ===");
    
    return 0;
}
