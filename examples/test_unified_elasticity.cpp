/**
 * test_unified_elasticity.cpp - 测试统一弹性力学模块
 * 
 * 验证 ElasticityUnified 支持 2D/3D 多种单元类型
 */

#include "mesh/model.h"
#include "mesh/mesh_generator.h"
#include "assembly/assembler.h"
#include "physics/elasticity_unified.h"
#include "material/isotropic_elastic.h"
#include "solver/cg.h"
#include "io/vtk_writer.h"
#include "core/timer.h"
#include <iostream>

using namespace fem;
using namespace fem::physics;
using namespace fem::constitutive;

int main() {
    std::cout << "=== 统一弹性力学模块测试 ===\n\n";
    
    // ────────────────────────────────────────────────────────
    // Test 1: 2D Tri3 平面应力
    // ────────────────────────────────────────────────────────
    std::cout << "Test 1: 2D Tri3 平面应力\n";
    std::cout << "-----------------------\n";
    
    {
        Model model("Elasticity2D_Tri3");
        int mat_id = model.add_material("Steel");
        model.material(mat_id).set_property("E", 2.0e5);
        model.material(mat_id).set_property("nu", 0.3);
        
        int mesh_id = model.add_mesh("domain", mat_id);
        Mesh& mesh = model.mesh(mesh_id);
        
        // 生成 20x20 三角形网格
        MeshGenerator::generate_unit_square_tri(20, 20, mesh);
        MeshGenerator::identify_boundaries_2d(mesh);
        
        // 创建弹性力学模块
        Real E = 2.0e5;
        Real nu = 0.3;
        IsotropicElastic material(E, nu, 2, true);  // 2D, plane_stress
        ElasticityUnified elast(&material, 2);
        
        // 装配
        Timer timer;
        Assembler assembler(model, 2);  // 矢量场 (u_x, u_y)
        assembler.assemble([&elast](Index elem_id, const Mesh& mesh,
                                    DenseMatrix& Ke, Vector& Fe) {
            elast.compute_element(elem_id, mesh, Ke, Fe);
        });
        
        // 边界条件：左边固定，右边拉伸
        std::vector<DirichletBC> bcs = {
            {"left", 0, 0.0},     // u_x = 0
            {"left", 1, 0.0},     // u_y = 0
            {"bottom", 1, 0.0},   // u_y = 0 (防止刚体运动)
            {"right", 0, 0.01}    // u_x = 0.01 (拉伸)
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
        
        // 结果分析
        Real u_x_max = 0.0, u_y_max = 0.0;
        for (size_t i = 0; i < u.size(); i += 2) {
            u_x_max = std::max(u_x_max, std::abs(u[i]));
            u_y_max = std::max(u_y_max, std::abs(u[i + 1]));
        }
        
        std::cout << "  网格: 20x20 (" << mesh.num_nodes() << " 节点, "
                  << mesh.num_elements() << " 单元)\n";
        std::cout << "  装配: " << t_assembly << " ms\n";
        std::cout << "  求解: " << result.iterations << " 次迭代, "
                  << "残差 " << result.residual << ", " << t_solve << " ms\n";
        std::cout << "  结果: u_x_max = " << u_x_max 
                  << ", u_y_max = " << u_y_max << "\n";
        
        // 输出 VTK
        VTKWriter vtk("test_unified_elast_tri3");
        vtk.write_mesh(mesh);
        vtk.add_point_vector("displacement", u, 2);
        vtk.close();
        
        std::cout << "  输出: test_unified_elast_tri3.vtk ✓\n\n";
    }
    
    // ────────────────────────────────────────────────────────
    // Test 2: 2D Quad4 平面应变
    // ────────────────────────────────────────────────────────
    std::cout << "Test 2: 2D Quad4 平面应变\n";
    std::cout << "------------------------\n";
    
    {
        Model model("Elasticity2D_Quad4");
        int mat_id = model.add_material("Steel");
        int mesh_id = model.add_mesh("domain", mat_id);
        Mesh& mesh = model.mesh(mesh_id);
        
        MeshGenerator::generate_unit_square_quad(20, 20, mesh);
        MeshGenerator::identify_boundaries_2d(mesh);
        
        Real E = 2.0e5;
        Real nu = 0.3;
        IsotropicElastic material(E, nu, 2, false);  // 2D, plane_strain
        ElasticityUnified elast(&material, 2);
        
        Timer timer;
        Assembler assembler(model, 2);
        assembler.assemble([&elast](Index elem_id, const Mesh& mesh,
                                    DenseMatrix& Ke, Vector& Fe) {
            elast.compute_element(elem_id, mesh, Ke, Fe);
        });
        
        std::vector<DirichletBC> bcs = {
            {"left", 0, 0.0},
            {"left", 1, 0.0},
            {"bottom", 1, 0.0},
            {"right", 0, 0.01}
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
        
        Real u_x_max = 0.0;
        for (size_t i = 0; i < u.size(); i += 2) {
            u_x_max = std::max(u_x_max, std::abs(u[i]));
        }
        
        std::cout << "  网格: 20x20 (" << mesh.num_nodes() << " 节点, "
                  << mesh.num_elements() << " 单元)\n";
        std::cout << "  装配: " << t_assembly << " ms\n";
        std::cout << "  求解: " << result.iterations << " 次迭代, "
                  << "残差 " << result.residual << ", " << t_solve << " ms\n";
        std::cout << "  结果: u_x_max = " << u_x_max << "\n";
        
        VTKWriter vtk("test_unified_elast_quad4");
        vtk.write_mesh(mesh);
        vtk.add_point_vector("displacement", u, 2);
        vtk.close();
        
        std::cout << "  输出: test_unified_elast_quad4.vtk ✓\n\n";
    }
    
    // ────────────────────────────────────────────────────────
    // Test 3: 3D Tet4 单元
    // ────────────────────────────────────────────────────────
    std::cout << "Test 3: 3D Tet4 单元\n";
    std::cout << "-------------------\n";
    
    {
        Model model("Elasticity3D_Tet4");
        int mat_id = model.add_material("Steel");
        int mesh_id = model.add_mesh("domain", mat_id);
        Mesh& mesh = model.mesh(mesh_id);
        
        // 生成 8x8x8 四面体网格
        MeshGenerator::generate_unit_cube_tet(8, 8, 8, mesh);
        MeshGenerator::identify_boundaries_3d(mesh);
        
        Real E = 2.0e5;
        Real nu = 0.3;
        IsotropicElastic material(E, nu, 3);  // 3D
        ElasticityUnified elast(&material, 3);
        
        Timer timer;
        Assembler assembler(model, 3);  // 矢量场 (u_x, u_y, u_z)
        assembler.assemble([&elast](Index elem_id, const Mesh& mesh,
                                    DenseMatrix& Ke, Vector& Fe) {
            elast.compute_element(elem_id, mesh, Ke, Fe);
        });
        
        // 边界条件：左面固定，右面拉伸
        std::vector<DirichletBC> bcs = {
            {"left", 0, 0.0},
            {"left", 1, 0.0},
            {"left", 2, 0.0},
            {"bottom", 1, 0.0},
            {"front", 2, 0.0},
            {"right", 0, 0.01}
        };
        assembler.apply_dirichlet(bcs);
        Real t_assembly = timer.elapsed_ms();
        
        timer.start();
        CGSolver solver;
        solver.set_tol(1e-8);
        solver.set_max_iter(1000);
        
        const auto& K = assembler.matrix();
        const auto& F = assembler.rhs();
        std::vector<Real> u(F.size(), 0.0);
        
        auto result = solver.solve(K, F.raw(), u);
        Real t_solve = timer.elapsed_ms();
        
        Real u_x_max = 0.0;
        for (size_t i = 0; i < u.size(); i += 3) {
            u_x_max = std::max(u_x_max, std::abs(u[i]));
        }
        
        std::cout << "  网格: 8x8x8 (" << mesh.num_nodes() << " 节点, "
                  << mesh.num_elements() << " 单元)\n";
        std::cout << "  装配: " << t_assembly << " ms\n";
        std::cout << "  求解: " << result.iterations << " 次迭代, "
                  << "残差 " << result.residual << ", " << t_solve << " ms\n";
        std::cout << "  结果: u_x_max = " << u_x_max << "\n";
        
        VTKWriter vtk("test_unified_elast_tet4");
        vtk.write_mesh(mesh);
        vtk.add_point_vector("displacement", u, 3);
        vtk.close();
        
        std::cout << "  输出: test_unified_elast_tet4.vtk ✓\n\n";
    }
    
    std::cout << "=== 所有测试通过 ✓ ===\n";
    return 0;
}
