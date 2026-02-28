/**
 * heat_2d.cpp - 2D 热传导示例 (使用新 physics::HeatConduction)
 * 
 * 问题: -∇·(k∇u) = Q
 * 
 * 边界条件: u = 0 (全边界 Dirichlet)
 * 导热系数: k = 1.0
 * 热源: Q = 10.0
 * 
 * 使用新架构: Model + Mesh + Assembler + HeatConduction
 */

#include "mesh/model.h"
#include "mesh/mesh_generator.h"
#include "assembly/assembler.h"
#include "physics/heat.h"
#include "math/cg.h"
#include "io/vtk_writer.h"
#include "core/timer.h"
#include "core/logger.h"

using namespace fem;
using namespace fem::physics;

int main() {
    FEM_INFO("=== 2D Heat Conduction Demo (new physics module) ===");
    
    // ── 1. 创建模型 ──
    Model model("Heat 2D");
    
    int mat_id = model.add_material("Thermal");
    model.material(mat_id).set_property("k", 1.0);
    model.material(mat_id).set_property("Q", 10.0);
    
    int mesh_id = model.add_mesh("domain", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    // 生成 30x30 网格
    MeshGenerator::generate_unit_square_tri(30, 30, mesh);
    MeshGenerator::identify_boundaries_2d(mesh);
    
    FEM_INFO("Mesh: " + std::to_string(mesh.num_nodes()) + " nodes, " + 
             std::to_string(mesh.num_elements()) + " elements");
    
    // ── 2. 创建物理模块 ──
    Timer timer;
    timer.start();
    
    Real k = mesh.material()->property("k", 1.0);
    Real Q = mesh.material()->property("Q", 10.0);
    
    HeatConduction heat(k, Q);
    
    // ── 3. 创建 Assembler 并装配 ──
    Assembler assembler(model, 1);  // 标量场
    
    // 使用 lambda 捕获 heat 物理模块
    auto elem_func = [&heat](Index elem_id, const Mesh& mesh,
                             DenseMatrix& Ke, Vector& Fe) {
        heat.compute_element(elem_id, mesh, Ke, Fe);
    };
    
    assembler.assemble(elem_func);
    
    FEM_INFO("Assembly done in: " + std::to_string(timer.elapsed_s()) + "s");
    
    // ── 4. 施加边界条件 ──
    timer.start();
    
    std::vector<DirichletBC> bcs;
    bcs.push_back({"left", 0, 0.0});
    bcs.push_back({"right", 0, 0.0});
    bcs.push_back({"bottom", 0, 0.0});
    bcs.push_back({"top", 0, 0.0});
    
    assembler.apply_dirichlet(bcs);
    
    FEM_INFO("Boundary conditions applied in: " + std::to_string(timer.elapsed_s()) + "s");
    
    // ── 5. 求解 ──
    timer.start();
    
    SparseMatrixCSR K = assembler.matrix();
    const Vector& F = assembler.rhs();
    
    std::vector<Real> F_std = F.raw();
    std::vector<Real> u_std(F.size(), 0.0);
    
    CGSolver solver;
    solver.set_tol(1e-8);
    solver.set_max_iter(2000);
    
    auto result = solver.solve(K, F_std, u_std);
    
    FEM_INFO("Solve completed in: " + std::to_string(timer.elapsed_s()) + "s");
    FEM_INFO("CG convergence: " + std::string(result.converged ? "YES" : "NO") + 
             ", residual=" + fmt_sci(result.residual));
    
    if (result.converged) {
        // ── 6. 输出结果 ──
        Real max_u = 0.0;
        std::size_t max_idx = 0;
        for (std::size_t i = 0; i < u_std.size(); ++i) {
            if (std::abs(u_std[i]) > max_u) {
                max_u = std::abs(u_std[i]);
                max_idx = i;
            }
        }
        
        const auto& max_coord = mesh.node(max_idx).coords();
        FEM_INFO("Max temperature: u = " + fmt_sci(max_u) + 
                 " at (" + fmt_sci(max_coord[0]) + ", " + fmt_sci(max_coord[1]) + ")");
        
        // 理论值: Q * L^2 / (8*k) = 10 * 1 / 8 = 1.25 (近似)
        FEM_INFO("Theoretical max (approx): 1.25");
        
        // ── 7. 输出 VTK 文件 ──
        FEM_INFO("\n=== Exporting VTK ===");
        try {
            VTKWriter vtk("heat_2d_result");
            vtk.write_mesh(mesh);
            vtk.add_point_scalar("temperature", u_std);
            vtk.close();
            FEM_INFO("VTK file written: heat_2d_result.vtk");
        } catch (const std::exception& e) {
            FEM_WARN("VTK export failed: " + std::string(e.what()));
        }
        
        FEM_INFO("Demo completed successfully!");
    } else {
        FEM_ERROR("Solver failed to converge");
        return 1;
    }
    
    return 0;
}
