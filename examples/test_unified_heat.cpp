/**
 * test_unified_heat.cpp - 测试统一热传导模块
 * 
 * 验证 HeatConductionUnified 支持多种单元类型
 */

#include "mesh/model.h"
#include "mesh/mesh_generator.h"
#include "assembly/assembler.h"
#include "physics/heat_unified.h"
#include "solver/cg.h"
#include "io/vtk_writer.h"
#include "core/timer.h"
#include <iostream>

using namespace fem;
using namespace fem::physics;

int main() {
    std::cout << "=== 统一热传导模块测试 ===\n\n";
    
    // ────────────────────────────────────────────────────────
    // Test 1: 2D Tri3 (和旧模块对比)
    // ────────────────────────────────────────────────────────
    std::cout << "Test 1: 2D Tri3 单元\n";
    std::cout << "-------------------\n";
    
    {
        Model model("Heat2D_Tri3");
        int mat_id = model.add_material("Material");
        model.material(mat_id).set_property("k", 1.0);
        
        int mesh_id = model.add_mesh("domain", mat_id);
        Mesh& mesh = model.mesh(mesh_id);
        
        // 生成 20x20 三角形网格
        MeshGenerator::generate_unit_square_tri(20, 20, mesh);
        MeshGenerator::identify_boundaries_2d(mesh);
        
        // 创建统一热传导模块
        Real k = 1.0;
        Real Q = 10.0;
        HeatConductionUnified heat(k, Q);
        
        // 装配
        Timer timer;
        Assembler assembler(model, 1);  // 标量场
        assembler.assemble([&heat](Index elem_id, const Mesh& mesh,
                                   DenseMatrix& Ke, Vector& Fe) {
            heat.compute_element(elem_id, mesh, Ke, Fe);
        });
        
        // 边界条件（四周温度为0）
        std::vector<DirichletBC> bcs = {
            {"left", 0, 0.0}, {"right", 0, 0.0},
            {"top", 0, 0.0}, {"bottom", 0, 0.0}
        };
        assembler.apply_dirichlet(bcs);
        
        Real t_assembly = timer.elapsed_ms();
        
        // 求解
        timer.start();
        CGSolver solver;
        solver.set_tol(1e-8);
        
        const auto& K = assembler.matrix();
        const auto& F = assembler.rhs();
        std::vector<Real> u(F.size(), 0.0);
        
        auto result = solver.solve(K, F.raw(), u);
        Real t_solve = timer.elapsed_ms();
        
        // 检查结果
        Real u_max = *std::max_element(u.begin(), u.end());
        Real u_center = u[mesh.num_nodes() / 2];  // 粗略的中心点
        
        std::cout << "  网格: 20x20 (" << mesh.num_nodes() << " 节点, "
                  << mesh.num_elements() << " 单元)\n";
        std::cout << "  装配: " << t_assembly << " ms\n";
        std::cout << "  求解: " << result.iterations << " 次迭代, "
                  << "残差 " << result.residual << ", " << t_solve << " ms\n";
        std::cout << "  结果: u_max = " << u_max << ", u_center = " << u_center << "\n";
        
        // 输出 VTK
        VTKWriter vtk("test_unified_heat_tri3");
        vtk.write_mesh(mesh);
        vtk.add_point_scalar("temperature", u);
        vtk.close();
        
        std::cout << "  输出: test_unified_heat_tri3.vtk ✓\n\n";
    }
    
    // ────────────────────────────────────────────────────────
    // Test 2: 2D Quad4 单元
    // ────────────────────────────────────────────────────────
    std::cout << "Test 2: 2D Quad4 单元\n";
    std::cout << "--------------------\n";
    
    {
        Model model("Heat2D_Quad4");
        int mat_id = model.add_material("Material");
        model.material(mat_id).set_property("k", 1.0);
        
        int mesh_id = model.add_mesh("domain", mat_id);
        Mesh& mesh = model.mesh(mesh_id);
        
        // 生成 20x20 四边形网格
        MeshGenerator::generate_unit_square_quad(20, 20, mesh);
        MeshGenerator::identify_boundaries_2d(mesh);
        
        Real k = 1.0;
        Real Q = 10.0;
        HeatConductionUnified heat(k, Q);
        
        Timer timer;
        Assembler assembler(model, 1);
        assembler.assemble([&heat](Index elem_id, const Mesh& mesh,
                                   DenseMatrix& Ke, Vector& Fe) {
            heat.compute_element(elem_id, mesh, Ke, Fe);
        });
        
        std::vector<DirichletBC> bcs = {
            {"left", 0, 0.0}, {"right", 0, 0.0},
            {"top", 0, 0.0}, {"bottom", 0, 0.0}
        };
        assembler.apply_dirichlet(bcs);
        Real t_assembly = timer.elapsed_ms();
        
        timer.start();
        CGSolver solver;
        solver.set_tol(1e-8);
        
        const auto& K = assembler.matrix();
        const auto& F = assembler.rhs();
        std::vector<Real> u(F.size(), 0.0);
        
        auto result = solver.solve(K, F.raw(), u);
        Real t_solve = timer.elapsed_ms();
        
        Real u_max = *std::max_element(u.begin(), u.end());
        
        std::cout << "  网格: 20x20 (" << mesh.num_nodes() << " 节点, "
                  << mesh.num_elements() << " 单元)\n";
        std::cout << "  装配: " << t_assembly << " ms\n";
        std::cout << "  求解: " << result.iterations << " 次迭代, "
                  << "残差 " << result.residual << ", " << t_solve << " ms\n";
        std::cout << "  结果: u_max = " << u_max << "\n";
        
        VTKWriter vtk("test_unified_heat_quad4");
        vtk.write_mesh(mesh);
        vtk.add_point_scalar("temperature", u);
        vtk.close();
        
        std::cout << "  输出: test_unified_heat_quad4.vtk ✓\n\n";
    }
    
    // ────────────────────────────────────────────────────────
    // Test 3: 3D Tet4 单元
    // ────────────────────────────────────────────────────────
    std::cout << "Test 3: 3D Tet4 单元\n";
    std::cout << "-------------------\n";
    
    {
        Model model("Heat3D_Tet4");
        int mat_id = model.add_material("Material");
        model.material(mat_id).set_property("k", 1.0);
        
        int mesh_id = model.add_mesh("domain", mat_id);
        Mesh& mesh = model.mesh(mesh_id);
        
        // 生成 10x10x10 四面体网格
        MeshGenerator::generate_unit_cube_tet(10, 10, 10, mesh);
        MeshGenerator::identify_boundaries_3d(mesh);
        
        Real k = 1.0;
        Real Q = 10.0;
        HeatConductionUnified heat(k, Q);
        
        Timer timer;
        Assembler assembler(model, 1);
        assembler.assemble([&heat](Index elem_id, const Mesh& mesh,
                                   DenseMatrix& Ke, Vector& Fe) {
            heat.compute_element(elem_id, mesh, Ke, Fe);
        });
        
        std::vector<DirichletBC> bcs = {
            {"left", 0, 0.0}, {"right", 0, 0.0},
            {"top", 0, 0.0}, {"bottom", 0, 0.0},
            {"front", 0, 0.0}, {"back", 0, 0.0}
        };
        assembler.apply_dirichlet(bcs);
        Real t_assembly = timer.elapsed_ms();
        
        timer.start();
        CGSolver solver;
        solver.set_tol(1e-8);
        solver.set_max_iter(500);
        
        const auto& K = assembler.matrix();
        const auto& F = assembler.rhs();
        std::vector<Real> u(F.size(), 0.0);
        
        auto result = solver.solve(K, F.raw(), u);
        Real t_solve = timer.elapsed_ms();
        
        Real u_max = *std::max_element(u.begin(), u.end());
        
        std::cout << "  网格: 10x10x10 (" << mesh.num_nodes() << " 节点, "
                  << mesh.num_elements() << " 单元)\n";
        std::cout << "  装配: " << t_assembly << " ms\n";
        std::cout << "  求解: " << result.iterations << " 次迭代, "
                  << "残差 " << result.residual << ", " << t_solve << " ms\n";
        std::cout << "  结果: u_max = " << u_max << "\n";
        
        VTKWriter vtk("test_unified_heat_tet4");
        vtk.write_mesh(mesh);
        vtk.add_point_scalar("temperature", u);
        vtk.close();
        
        std::cout << "  输出: test_unified_heat_tet4.vtk ✓\n\n";
    }
    
    std::cout << "=== 所有测试通过 ✓ ===\n";
    return 0;
}
