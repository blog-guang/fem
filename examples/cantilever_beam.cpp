/**
 * cantilever_beam.cpp - 悬臂梁示例 (Neumann 边界条件)
 * 
 * 问题: 2D 悬臂梁在顶部受均布载荷
 * 
 * 边界条件:
 * - 左边界固定: u_x = u_y = 0 (Dirichlet)
 * - 顶部载荷: t_y = -p (向下, Neumann)
 * - 其他边界自由
 * 
 * 材料参数:
 * - E = 1e6 (杨氏模量)
 * - ν = 0.3 (泊松比)
 * 
 * 几何: 4x1 矩形梁 (长宽比 4:1)
 */

#include "mesh/model.h"
#include "mesh/mesh_generator.h"
#include "assembly/assembler.h"
#include "physics/elasticity_unified.h"
#include "material/isotropic_elastic.h"
#include "math/cg.h"
#include "io/vtk_writer.h"
#include "core/timer.h"
#include "core/logger.h"

using namespace fem;
using namespace fem::physics;

int main() {
    FEM_INFO("=== Cantilever Beam with Neumann BC ===");
    
    // ── 1. 创建模型 ──
    Model model("Cantilever");
    
    int mat_id = model.add_material("Steel");
    model.material(mat_id).set_property("E", 1e6);   // 杨氏模量
    model.material(mat_id).set_property("nu", 0.3);  // 泊松比
    
    int mesh_id = model.add_mesh("beam", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    // 生成 40x10 网格 (长4 x 高1，细化到保证精度)
    MeshGenerator::generate_unit_square_tri(40, 10, mesh);
    
    // 缩放到 4x1 矩形
    for (std::size_t i = 0; i < mesh.num_nodes(); ++i) {
        auto& coords = mesh.node(i).coords();
        coords[0] *= 4.0;  // 长度 4
        coords[1] *= 1.0;  // 高度 1
    }
    
    MeshGenerator::identify_boundaries_2d(mesh);
    
    FEM_INFO("Mesh: " + std::to_string(mesh.num_nodes()) + " nodes, " + 
             std::to_string(mesh.num_elements()) + " elements");
    
    // ── 2. 创建物理模块 ──
    Timer timer;
    timer.start();
    
    Real E = mesh.material()->property("E", 1e6);
    Real nu = mesh.material()->property("nu", 0.3);
    
    // 使用线弹性材料
    constitutive::IsotropicElastic material(E, nu, 2);
    ElasticityUnified elast(&material, 2);
    
    // ── 3. 创建 Assembler 并装配 ──
    Assembler assembler(model, 2);  // 矢量场 (u_x, u_y)
    
    auto elem_func = [&elast](Index elem_id, const Mesh& mesh,
                              DenseMatrix& Ke, Vector& Fe) {
        elast.compute_element(elem_id, mesh, Ke, Fe);
    };
    
    assembler.assemble(elem_func);
    
    timer.stop();
    FEM_INFO("Assembly: " + std::to_string(timer.elapsed_ms()) + "ms");
    
    // ── 4. 施加 Dirichlet 边界条件 (左边固定) ──
    std::vector<DirichletBC> dirichlet_bcs = {
        {"left", 0, 0.0},   // u_x = 0
        {"left", 1, 0.0}    // u_y = 0
    };
    
    assembler.apply_dirichlet(dirichlet_bcs);
    
    // ── 5. 施加 Neumann 边界条件 (顶部均布载荷) ──
    Real pressure = -10.0;  // 向下的压力
    
    std::vector<NeumannBC> neumann_bcs = {
        {"top", 1, pressure}  // t_y = -10 (y方向, 向下)
    };
    
    assembler.apply_neumann(neumann_bcs);
    
    // ── 6. 求解 ──
    timer.start();
    
    const auto& K = assembler.matrix();
    const auto& F = assembler.rhs();
    
    Vector u(F.size(), 0.0);
    
    CGSolver solver;
    solver.set_tol(1e-8);
    solver.set_max_iter(10000);
    
    auto result = solver.solve(K, F, u);
    
    timer.stop();
    
    if (!result.converged) {
        FEM_ERROR("Solver failed to converge");
        return 1;
    }
    
    FEM_INFO("Solve: " + std::to_string(result.iterations) + " iterations, " +
             std::to_string(timer.elapsed_ms()) + "ms");
    
    // ── 7. 后处理 ──
    
    // 找最大位移
    Real max_u_y = 0.0;
    std::size_t max_idx = 0;
    for (std::size_t i = 0; i < mesh.num_nodes(); ++i) {
        Real u_y = u[i * 2 + 1];
        if (std::abs(u_y) > std::abs(max_u_y)) {
            max_u_y = u_y;
            max_idx = i;
        }
    }
    
    const auto& max_coord = mesh.node(max_idx).coords();
    
    FEM_INFO("Max displacement: |u_y| = " + fmt_sci(std::abs(max_u_y)) + 
             " at (" + fmt_sci(max_coord[0]) + ", " + fmt_sci(max_coord[1]) + ")");
    
    // 理论解 (悬臂梁自由端最大挠度)
    // δ_max = (p * L^4) / (8 * E * I)
    // I = b*h³/12 (矩形梁惯性矩)
    Real L = 4.0;
    Real h = 1.0;
    Real b = 1.0;  // 单位厚度 (平面应力)
    Real I = b * h * h * h / 12.0;
    Real delta_theory = std::abs(pressure) * L * L * L * L / (8.0 * E * I);
    
    FEM_INFO("Theoretical: " + fmt_sci(delta_theory));
    FEM_INFO("Error: " + fmt_sci(std::abs(std::abs(max_u_y) - delta_theory) / delta_theory * 100.0) + "%");
    
    // ── 8. 输出 VTK 文件 ──
    VTKWriter vtk("cantilever_beam_result");
    vtk.write_mesh(mesh);
    vtk.add_point_vector("displacement", u, 2);
    vtk.close();
    FEM_INFO("VTK: cantilever_beam_result.vtk");
    
    FEM_INFO("=== Demo completed successfully! ===");
    
    return 0;
}
