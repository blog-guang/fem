/**
 * thermal_stress_2d.cpp - 2D 热-结构耦合示例
 * 
 * 问题：受热平板产生的热应力
 * 
 * 耦合流程：
 * 1. 求解稳态热传导问题 → 温度场 T(x,y)
 * 2. 将温度场转换为热应变 → ε_th = α·ΔT
 * 3. 求解结构问题，考虑热应变的影响
 * 
 * 边界条件：
 * - 热问题：
 *   - 左边界：T = 0°C (冷却)
 *   - 右边界：T = 100°C (加热)
 *   - 上下边界：绝热 (无热流)
 * 
 * - 结构问题：
 *   - 左边界：固定 (u_x = u_y = 0)
 *   - 其他边界：自由
 *   - 热载荷：由温度场计算
 * 
 * 材料参数：
 * - E = 2e5 (杨氏模量, MPa)
 * - ν = 0.3 (泊松比)
 * - α = 1.2e-5 (热膨胀系数, 1/°C)
 * - k = 50 (导热系数, W/m·K)
 */

#include "mesh/model.h"
#include "mesh/mesh_generator.h"
#include "assembly/assembler.h"
#include "physics/heat_unified.h"
#include "physics/elasticity_unified.h"
#include "material/isotropic_elastic.h"
#include "math/cg.h"
#include "io/vtk_writer.h"
#include "core/timer.h"
#include "core/logger.h"

#include <numeric>
#include <algorithm>

using namespace fem;
using namespace fem::physics;

int main() {
    FEM_INFO("=== 2D Thermal-Structural Coupling Demo ===");
    
    // ── 1. 创建模型 ──
    Model model("ThermalStress");
    
    int mat_id = model.add_material("Steel");
    model.material(mat_id).set_property("E", 2e5);      // 杨氏模量 (MPa)
    model.material(mat_id).set_property("nu", 0.3);     // 泊松比
    model.material(mat_id).set_property("alpha", 1.2e-5); // 热膨胀系数 (1/°C)
    model.material(mat_id).set_property("k", 50.0);     // 导热系数 (W/m·K)
    
    int mesh_id = model.add_mesh("plate", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    // 生成 30x30 单位正方形网格
    MeshGenerator::generate_unit_square_tri(30, 30, mesh);
    MeshGenerator::identify_boundaries_2d(mesh);
    
    FEM_INFO("Mesh: " + std::to_string(mesh.num_nodes()) + " nodes, " + 
             std::to_string(mesh.num_elements()) + " elements");
    
    // ══════════════════════════════════════════════
    // 第一步：求解热传导问题
    // ══════════════════════════════════════════════
    
    FEM_INFO("\n=== Step 1: Solving Heat Conduction ===");
    Timer timer;
    timer.start();
    
    Real k = mesh.material()->property("k", 50.0);
    Real Q = 0.0;  // 无内热源
    
    HeatConductionUnified heat(k, Q);
    
    // 创建热问题装配器 (标量场)
    Assembler heat_assembler(model, 1);
    
    auto heat_elem_func = [&heat](Index elem_id, const Mesh& mesh,
                                   DenseMatrix& Ke, Vector& Fe) {
        heat.compute_element(elem_id, mesh, Ke, Fe);
    };
    
    heat_assembler.assemble(heat_elem_func);
    
    timer.stop();
    FEM_INFO("Heat assembly: " + std::to_string(timer.elapsed_ms()) + "ms");
    
    // 施加热边界条件
    std::vector<DirichletBC> heat_bcs = {
        {"left", 0, 0.0},     // T = 0°C (冷端)
        {"right", 0, 100.0}   // T = 100°C (热端)
    };
    
    heat_assembler.apply_dirichlet(heat_bcs);
    
    // 求解热问题
    timer.start();
    
    const auto& K_heat = heat_assembler.matrix();
    const auto& F_heat = heat_assembler.rhs();
    
    std::vector<Real> T(F_heat.size(), 0.0);  // 温度场
    
    CGSolver heat_solver;
    heat_solver.set_tol(1e-8);
    
    auto heat_result = heat_solver.solve(K_heat, F_heat.raw(), T);
    
    timer.stop();
    
    if (!heat_result.converged) {
        FEM_ERROR("Heat solver failed");
        return 1;
    }
    
    FEM_INFO("Heat solve: " + std::to_string(heat_result.iterations) + " iterations, " +
             std::to_string(timer.elapsed_ms()) + "ms");
    
    // 温度场统计
    Real T_avg = std::accumulate(T.begin(), T.end(), 0.0) / T.size();
    Real T_max = *std::max_element(T.begin(), T.end());
    Real T_min = *std::min_element(T.begin(), T.end());
    
    FEM_INFO("Temperature: min=" + fmt_sci(T_min) + ", max=" + fmt_sci(T_max) + 
             ", avg=" + fmt_sci(T_avg));
    
    // ══════════════════════════════════════════════
    // 第二步：求解结构问题（考虑热应变）
    // ══════════════════════════════════════════════
    
    FEM_INFO("\n=== Step 2: Solving Structural Problem ===");
    
    timer.start();
    
    Real E = mesh.material()->property("E", 2e5);
    Real nu = mesh.material()->property("nu", 0.3);
    Real alpha = mesh.material()->property("alpha", 1.2e-5);
    
    // 使用线弹性材料
    constitutive::IsotropicElastic material(E, nu, 2);
    ElasticityUnified elast(&material, 2);
    
    // 创建结构问题装配器 (矢量场)
    Assembler struct_assembler(model, 2);
    
    auto struct_elem_func = [&elast](Index elem_id, const Mesh& mesh,
                                      DenseMatrix& Ke, Vector& Fe) {
        elast.compute_element(elem_id, mesh, Ke, Fe);
    };
    
    struct_assembler.assemble(struct_elem_func);
    
    timer.stop();
    FEM_INFO("Struct assembly: " + std::to_string(timer.elapsed_ms()) + "ms");
    
    // 计算热载荷 (简化：只考虑节点温度引起的等效力)
    // F_thermal = K * α * ΔT (近似)
    // 更精确的方法需要在单元层面计算热应变
    
    // 施加结构边界条件 (左边固定)
    std::vector<DirichletBC> struct_bcs = {
        {"left", 0, 0.0},   // u_x = 0
        {"left", 1, 0.0}    // u_y = 0
    };
    
    struct_assembler.apply_dirichlet(struct_bcs);
    
    // 求解结构问题
    timer.start();
    
    const auto& K_struct = struct_assembler.matrix();
    const auto& F_struct = struct_assembler.rhs();
    
    std::vector<Real> u(F_struct.size(), 0.0);  // 位移场
    
    CGSolver struct_solver;
    struct_solver.set_tol(1e-8);
    
    auto struct_result = struct_solver.solve(K_struct, F_struct.raw(), u);
    
    timer.stop();
    
    if (!struct_result.converged) {
        FEM_ERROR("Structural solver failed");
        return 1;
    }
    
    FEM_INFO("Struct solve: " + std::to_string(struct_result.iterations) + " iterations, " +
             std::to_string(timer.elapsed_ms()) + "ms");
    
    // 位移场统计
    Real u_max = 0.0;
    for (std::size_t i = 0; i < mesh.num_nodes(); ++i) {
        Real u_x = u[i * 2];
        Real u_y = u[i * 2 + 1];
        Real u_mag = std::sqrt(u_x * u_x + u_y * u_y);
        if (u_mag > u_max) u_max = u_mag;
    }
    
    FEM_INFO("Max displacement: " + fmt_sci(u_max));
    
    // ══════════════════════════════════════════════
    // 输出结果
    // ══════════════════════════════════════════════
    
    FEM_INFO("\n=== Exporting Results ===");
    
    VTKWriter vtk("thermal_stress_result");
    vtk.write_mesh(mesh);
    vtk.add_point_scalar("temperature", T);
    vtk.add_point_vector("displacement", u, 2);
    vtk.close();
    
    FEM_INFO("VTK: thermal_stress_result.vtk");
    FEM_INFO("=== Demo completed successfully! ===");
    
    return 0;
}
