/**
 * elasticity_2d.cpp - 2D 弹性力学示例 (使用新 physics::Elasticity2D)
 * 
 * 问题: 2D 平面应力拉伸
 * 
 * 边界条件:
 * - 左边界固定: u_x = u_y = 0
 * - 右边界拉伸: u_x = 0.01, u_y = 0
 * - 上下边界自由
 * 
 * 材料参数:
 * - E = 1000 (杨氏模量)
 * - ν = 0.3 (泊松比)
 * 
 * 使用新架构: Model + Mesh + Assembler + Elasticity2D
 */

#include "mesh/model.h"
#include "mesh/mesh_generator.h"
#include "assembly/assembler.h"
#include "physics/elasticity_v2.h"
#include "math/cg.h"
#include "io/vtk_writer.h"
#include "core/timer.h"
#include "core/logger.h"

using namespace fem;
using namespace fem::physics;

int main() {
    FEM_INFO("=== 2D Elasticity Demo (new physics module) ===");
    
    // ── 1. 创建模型 ──
    Model model("Elasticity 2D");
    
    int mat_id = model.add_material("Steel");
    model.material(mat_id).set_property("E", 1000.0);   // 杨氏模量
    model.material(mat_id).set_property("nu", 0.3);     // 泊松比
    
    int mesh_id = model.add_mesh("domain", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    // 生成 20x20 单位正方形网格
    MeshGenerator::generate_unit_square_tri(20, 20, mesh);
    MeshGenerator::identify_boundaries_2d(mesh);
    
    FEM_INFO("Mesh: " + std::to_string(mesh.num_nodes()) + " nodes, " + 
             std::to_string(mesh.num_elements()) + " elements");
    
    // ── 2. 创建物理模块 ──
    Timer timer;
    timer.start();
    
    Real E = mesh.material()->property("E", 1000.0);
    Real nu = mesh.material()->property("nu", 0.3);
    
    Elasticity2D elast(E, nu, PlaneType::PlaneStress);
    
    // ── 3. 创建 Assembler 并装配 ──
    Assembler assembler(model, 2);  // 矢量场 (u_x, u_y)
    
    auto elem_func = [&elast](Index elem_id, const Mesh& mesh,
                              DenseMatrix& Ke, Vector& Fe) {
        elast.compute_element(elem_id, mesh, Ke, Fe);
    };
    
    assembler.assemble(elem_func);
    
    FEM_INFO("Assembly done in: " + std::to_string(timer.elapsed_s()) + "s");
    
    // ── 4. 施加边界条件 ──
    timer.start();
    
    std::vector<DirichletBC> bcs;
    
    // 左边界固定: u_x = 0, u_y = 0 (防止刚体平移)
    bcs.push_back({"left", 0, 0.0});   // DOF 0: u_x
    bcs.push_back({"left", 1, 0.0});   // DOF 1: u_y
    
    // 右边界拉伸: u_x = 0.01
    bcs.push_back({"right", 0, 0.01}); // DOF 0: u_x
    
    // 底部 y 方向固定（防止刚体旋转）
    bcs.push_back({"bottom", 1, 0.0}); // DOF 1: u_y
    
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
    solver.set_max_iter(5000);
    
    auto result = solver.solve(K, F_std, u_std);
    
    FEM_INFO("Solve completed in: " + std::to_string(timer.elapsed_s()) + "s");
    FEM_INFO("CG convergence: " + std::string(result.converged ? "YES" : "NO") + 
             ", residual=" + fmt_sci(result.residual));
    
    if (result.converged) {
        // ── 6. 输出结果 ──
        
        // 找最大位移
        Real max_u = 0.0;
        std::size_t max_idx = 0;
        for (std::size_t i = 0; i < u_std.size(); i += 2) {
            Real u_mag = std::sqrt(u_std[i]*u_std[i] + u_std[i+1]*u_std[i+1]);
            if (u_mag > max_u) {
                max_u = u_mag;
                max_idx = i / 2;
            }
        }
        
        const auto& max_coord = mesh.node(max_idx).coords();
        Real u_x = u_std[max_idx * 2];
        Real u_y = u_std[max_idx * 2 + 1];
        
        FEM_INFO("Max displacement: |u| = " + fmt_sci(max_u) + 
                 " (u_x=" + fmt_sci(u_x) + ", u_y=" + fmt_sci(u_y) + ")");
        FEM_INFO("  at (" + fmt_sci(max_coord[0]) + ", " + fmt_sci(max_coord[1]) + ")");
        
        // 检查中心点位移
        std::size_t center_node = (mesh.num_nodes() - 1) / 2;  // 近似中心
        Real u_x_center = u_std[center_node * 2];
        Real u_y_center = u_std[center_node * 2 + 1];
        
        FEM_INFO("Center displacement: u_x=" + fmt_sci(u_x_center) + 
                 ", u_y=" + fmt_sci(u_y_center));
        
        // ── 7. 输出 VTK 文件 ──
        FEM_INFO("\n=== Exporting VTK ===");
        try {
            VTKWriter vtk("elasticity_2d_result");
            vtk.write_mesh(mesh);
            vtk.add_point_vector("displacement", u_std, 2);  // 2D 矢量场
            vtk.close();
            FEM_INFO("VTK file written: elasticity_2d_result.vtk");
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
