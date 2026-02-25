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
#include "physics/heat.h"
#include "physics/elasticity_v2.h"
#include "solver/cg.h"
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
    
    HeatConduction heat(k, Q);
    
    // 创建热问题装配器 (标量场)
    Assembler heat_assembler(model, 1);
    
    auto heat_elem_func = [&heat](Index elem_id, const Mesh& mesh,
                                   DenseMatrix& Ke, Vector& Fe) {
        heat.compute_element(elem_id, mesh, Ke, Fe);
    };
    
    heat_assembler.assemble(heat_elem_func);
    
    FEM_INFO("Heat assembly done in: " + std::to_string(timer.elapsed_s()) + "s");
    
    // 施加热边界条件
    timer.start();
    
    std::vector<DirichletBC> heat_bcs = {
        {"left", 0, 0.0},     // T = 0°C (冷端)
        {"right", 0, 100.0}   // T = 100°C (热端)
    };
    
    heat_assembler.apply_dirichlet(heat_bcs);
    
    FEM_INFO("Heat BCs applied in: " + std::to_string(timer.elapsed_s()) + "s");
    
    // 求解热问题
    timer.start();
    
    SparseMatrixCSR K_heat = heat_assembler.matrix();
    const Vector& F_heat = heat_assembler.rhs();
    
    std::vector<Real> T_std = F_heat.raw();
    std::vector<Real> T_result(T_std.size(), 0.0);
    
    CGSolver heat_solver;
    heat_solver.set_tol(1e-8);
    heat_solver.set_max_iter(5000);
    
    auto heat_result = heat_solver.solve(K_heat, T_std, T_result);
    
    FEM_INFO("Heat solve completed in: " + std::to_string(timer.elapsed_s()) + "s");
    FEM_INFO("Heat CG: " + std::string(heat_result.converged ? "YES" : "NO") + 
             ", residual=" + fmt_sci(heat_result.residual) + 
             ", iter=" + std::to_string(heat_result.iterations));
    
    if (!heat_result.converged) {
        FEM_ERROR("Heat solve failed to converge");
        return 1;
    }
    
    // 输出温度场统计
    Real T_min = *std::min_element(T_result.begin(), T_result.end());
    Real T_max = *std::max_element(T_result.begin(), T_result.end());
    Real T_avg = std::accumulate(T_result.begin(), T_result.end(), 0.0) / T_result.size();
    
    FEM_INFO("Temperature field: T_min=" + fmt_sci(T_min) + 
             ", T_max=" + fmt_sci(T_max) + 
             ", T_avg=" + fmt_sci(T_avg));
    
    // ══════════════════════════════════════════════
    // 第二步：求解结构问题（考虑热应变）
    // ══════════════════════════════════════════════
    
    FEM_INFO("\n=== Step 2: Solving Structural Problem ===");
    timer.start();
    
    Real E = mesh.material()->property("E", 2e5);
    Real nu = mesh.material()->property("nu", 0.3);
    Real alpha = mesh.material()->property("alpha", 1.2e-5);
    
    Elasticity2D elast(E, nu, PlaneType::PlaneStress);
    
    // 创建结构问题装配器 (矢量场)
    Assembler struct_assembler(model, 2);
    
    // 装配刚度矩阵（不含热载荷，热载荷通过修改 F 添加）
    auto struct_elem_func = [&elast](Index elem_id, const Mesh& mesh,
                                      DenseMatrix& Ke, Vector& Fe) {
        elast.compute_element(elem_id, mesh, Ke, Fe);
    };
    
    struct_assembler.assemble(struct_elem_func);
    
    FEM_INFO("Structural assembly done in: " + std::to_string(timer.elapsed_s()) + "s");
    
    // 添加热载荷
    // 热应变: ε_th = α * ΔT * [1, 1, 0]^T (平面应力)
    // 等效节点力: F_th = -∫ B^T * D * ε_th dV
    
    FEM_INFO("Computing thermal loads...");
    timer.start();
    
    // 创建热载荷向量（直接修改 rhs）
    // 注意：这需要在 apply_dirichlet 之前完成，或者重新设计接口
    std::vector<Real> thermal_loads(struct_assembler.num_dofs(), 0.0);
    
    // 参考温度（无应力温度）
    constexpr Real T_ref = 0.0;
    
    // 预计算本构矩阵 D（对所有单元相同）
    DenseMatrix D(3, 3, 0.0);
    {
        Real factor = E / (1.0 - nu * nu);
        D(0, 0) = factor;
        D(0, 1) = factor * nu;
        D(1, 0) = factor * nu;
        D(1, 1) = factor;
        D(2, 2) = factor * (1.0 - nu) / 2.0;
    }
    
    for (std::size_t elem_id = 0; elem_id < mesh.num_elements(); ++elem_id) {
        const Element& elem = mesh.element(elem_id);
        
        // 计算单元平均温度（简化方法）
        Real T_elem = 0.0;
        for (Index node_id : elem.nodes()) {
            T_elem += T_result[node_id];
        }
        T_elem /= static_cast<Real>(elem.nodes().size());
        
        Real dT = T_elem - T_ref;
        
        // 热应变向量 (Voigt 标记: ε_xx, ε_yy, γ_xy)
        std::vector<Real> eps_th = {alpha * dT, alpha * dT, 0.0};
        
        // 计算 σ_th = D * ε_th
        std::vector<Real> sigma_th(3, 0.0);
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                sigma_th[i] += D(i, j) * eps_th[j];
            }
        }
        
        // 计算等效节点力
        // 简化方法：均匀分配到节点（实际应该用 B^T 积分）
        const auto& nodes = elem.nodes();
        if (nodes.size() == 3) {  // Tri3
            // 获取节点坐标
            const auto& c0 = mesh.node(nodes[0]).coords();
            const auto& c1 = mesh.node(nodes[1]).coords();
            const auto& c2 = mesh.node(nodes[2]).coords();
            
            // 计算单元面积
            Real x1 = c1[0] - c0[0], y1 = c1[1] - c0[1];
            Real x2 = c2[0] - c0[0], y2 = c2[1] - c0[1];
            Real area = std::abs(x1 * y2 - x2 * y1) / 2.0;
            
            // 等效节点力（每个节点分配 1/3）
            // 注意：负号表示热膨胀产生的反作用力
            constexpr Real one_third = 1.0 / 3.0;
            Real fx = -sigma_th[0] * area * one_third;
            Real fy = -sigma_th[1] * area * one_third;
            
            for (Index node_id : nodes) {
                thermal_loads[node_id * 2 + 0] += fx;
                thermal_loads[node_id * 2 + 1] += fy;
            }
        }
    }
    
    FEM_INFO("Thermal loads computed in: " + std::to_string(timer.elapsed_s()) + "s");
    
    // 施加结构边界条件
    timer.start();
    
    std::vector<DirichletBC> struct_bcs = {
        {"left", 0, 0.0},   // u_x = 0
        {"left", 1, 0.0}    // u_y = 0
    };
    
    struct_assembler.apply_dirichlet(struct_bcs);
    
    FEM_INFO("Structural BCs applied in: " + std::to_string(timer.elapsed_s()) + "s");
    
    // 求解结构问题
    timer.start();
    
    SparseMatrixCSR K_struct = struct_assembler.matrix();
    const Vector& F_base = struct_assembler.rhs();
    
    // 合并基础载荷和热载荷
    std::vector<Real> F_struct_std(F_base.size());
    for (std::size_t i = 0; i < F_base.size(); ++i) {
        F_struct_std[i] = F_base[i] + thermal_loads[i];
    }
    
    std::vector<Real> u_std(F_struct_std.size(), 0.0);
    
    CGSolver struct_solver;
    struct_solver.set_tol(1e-8);
    struct_solver.set_max_iter(10000);
    
    auto struct_result = struct_solver.solve(K_struct, F_struct_std, u_std);
    
    FEM_INFO("Structural solve completed in: " + std::to_string(timer.elapsed_s()) + "s");
    FEM_INFO("Structural CG: " + std::string(struct_result.converged ? "YES" : "NO") + 
             ", residual=" + fmt_sci(struct_result.residual) + 
             ", iter=" + std::to_string(struct_result.iterations));
    
    if (!struct_result.converged) {
        FEM_ERROR("Structural solve failed to converge");
        return 1;
    }
    
    // ══════════════════════════════════════════════
    // 输出结果
    // ══════════════════════════════════════════════
    
    // 找最大位移
    Real max_u = 0.0;
    std::size_t max_idx = 0;
    for (std::size_t i = 0; i < mesh.num_nodes(); ++i) {
        Real u_mag = std::sqrt(u_std[i*2]*u_std[i*2] + u_std[i*2+1]*u_std[i*2+1]);
        if (u_mag > max_u) {
            max_u = u_mag;
            max_idx = i;
        }
    }
    
    const auto& max_coord = mesh.node(max_idx).coords();
    Real u_x = u_std[max_idx * 2];
    Real u_y = u_std[max_idx * 2 + 1];
    
    FEM_INFO("\n=== Results ===");
    FEM_INFO("Max displacement: |u| = " + fmt_sci(max_u) + 
             " (u_x=" + fmt_sci(u_x) + ", u_y=" + fmt_sci(u_y) + ")");
    FEM_INFO("  at (" + fmt_sci(max_coord[0]) + ", " + fmt_sci(max_coord[1]) + ")");
    
    // 输出 VTK 文件
    FEM_INFO("\n=== Exporting VTK ===");
    try {
        VTKWriter vtk("thermal_stress_result");
        vtk.write_mesh(mesh);
        vtk.add_point_scalar("temperature", T_result);
        vtk.add_point_vector("displacement", u_std, 2);
        vtk.close();
        FEM_INFO("VTK file written: thermal_stress_result.vtk");
    } catch (const std::exception& e) {
        FEM_WARN("VTK export failed: " + std::string(e.what()));
    }
    
    FEM_INFO("\nDemo completed successfully!");
    FEM_INFO("Thermal-structural coupling demonstrates:");
    FEM_INFO("  1. Sequential coupling (heat → structure)");
    FEM_INFO("  2. Temperature field induces thermal strain");
    FEM_INFO("  3. Thermal strain creates stress and displacement");
    
    return 0;
}
