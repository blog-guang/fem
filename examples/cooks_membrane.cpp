/**
 * cooks_membrane.cpp - Cook's Membrane 经典几何非线性问题
 * 
 * 问题描述：
 * - 梯形板，受端部剪切载荷
 * - 平面应变条件
 * - 几何非线性大变形
 * 
 * 几何参数：
 * - 左边长度：44 mm
 * - 右边长度：16 mm
 * - 宽度：48 mm
 * - 厚度：1 mm（平面应变）
 * 
 * 材料参数：
 * - 杨氏模量：E = 250 Pa（软材料，突出非线性）
 * - 泊松比：ν = 0.4999（近不可压）
 * 
 * 边界条件：
 * - 左边固支（u = v = 0）
 * - 右边受剪切载荷（τ = 100 Pa，向上）
 * 
 * 参考解（文献值）：
 * - 右上角竖向位移：~23.96 mm（线性）
 * - 右上角竖向位移：~24.05 mm（几何非线性，P.Wriggers 2008）
 * 
 * ANSYS 验证数据（如有）：
 * - PLANE182, 8x8 网格
 * - NLGEOM,ON
 */

#include "mesh/model.h"
#include "mesh/mesh_generator.h"
#include "material/neo_hookean.h"
#include "physics/elasticity_nonlinear.h"
#include "solver/newton_raphson_solver.h"
#include "assembly/assembler.h"
#include "io/vtk_writer.h"
#include "core/logger.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace fem;
using namespace fem::constitutive;
using namespace fem::physics;
using namespace fem::solver;

/**
 * 生成 Cook's Membrane 网格
 * 
 * @param nx X 方向单元数
 * @param ny Y 方向单元数
 * @param mesh 目标网格
 */
void generate_cooks_membrane_mesh(int nx, int ny, Mesh& mesh) {
    // 几何参数
    Real width = 48.0;      // 宽度 mm
    Real h_left = 44.0;     // 左边高度 mm
    Real h_right = 16.0;    // 右边高度 mm
    
    // 生成节点
    std::vector<Index> node_map((nx + 1) * (ny + 1));
    
    for (int j = 0; j <= ny; j++) {
        for (int i = 0; i <= nx; i++) {
            Real xi = static_cast<Real>(i) / nx;
            Real eta = static_cast<Real>(j) / ny;
            
            // X 坐标：线性分布
            Real x = xi * width;
            
            // Y 坐标：梯形插值
            Real h_local = h_left * (1.0 - xi) + h_right * xi;
            Real y = eta * h_local;
            
            Vec3 coords = {x, y, 0.0};
            Index node_id = mesh.add_node(coords);
            node_map[j * (nx + 1) + i] = node_id;
        }
    }
    
    // 生成 Quad4 单元
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            Index n0 = node_map[j * (nx + 1) + i];
            Index n1 = node_map[j * (nx + 1) + (i + 1)];
            Index n2 = node_map[(j + 1) * (nx + 1) + (i + 1)];
            Index n3 = node_map[(j + 1) * (nx + 1) + i];
            
            mesh.add_element(ElementType::Quad4, {n0, n1, n2, n3});
        }
    }
    
    LOG_INFO("Generated Cook's Membrane mesh: " << nx << "x" << ny);
    LOG_INFO("  Nodes: " << mesh.num_nodes());
    LOG_INFO("  Elements: " << mesh.num_elements());
}

/**
 * 应用边界条件
 * 
 * @param mesh 网格
 * @param assembler 装配器
 * @param shear_stress 剪切应力
 * @return 外力向量
 */
Vector apply_boundary_conditions(
    const Mesh& mesh,
    Assembler& assembler,
    Real shear_stress)
{
    int dofs_per_node = 2;  // 平面问题
    Vector F_ext(mesh.num_nodes() * dofs_per_node, 0.0);
    
    // 左边固支（x = 0）
    for (Index node_id = 0; node_id < mesh.num_nodes(); node_id++) {
        const auto& coords = mesh.node(node_id).coords();
        
        if (std::abs(coords[0]) < 1e-6) {  // x ≈ 0
            assembler.apply_dirichlet(node_id, 0, 0.0);  // u = 0
            assembler.apply_dirichlet(node_id, 1, 0.0);  // v = 0
        }
    }
    
    // 右边剪切载荷（x = 48）
    // 分布载荷：F = τ * h * thickness
    // 假设厚度 = 1 mm
    Real thickness = 1.0;
    Real width = 48.0;
    
    // 找到右边节点（x ≈ 48）
    std::vector<Index> right_nodes;
    for (Index node_id = 0; node_id < mesh.num_nodes(); node_id++) {
        const auto& coords = mesh.node(node_id).coords();
        
        if (std::abs(coords[0] - width) < 1e-6) {
            right_nodes.push_back(node_id);
        }
    }
    
    // 按 y 坐标排序
    std::sort(right_nodes.begin(), right_nodes.end(),
        [&mesh](Index a, Index b) {
            return mesh.node(a).coords()[1] < mesh.node(b).coords()[1];
        });
    
    // 分配载荷（梯形积分）
    for (std::size_t i = 0; i < right_nodes.size(); i++) {
        Index node_id = right_nodes[i];
        Real dy = 0.0;
        
        if (i == 0 || i == right_nodes.size() - 1) {
            // 端点：半个段
            if (i == 0 && right_nodes.size() > 1) {
                Real y0 = mesh.node(right_nodes[i]).coords()[1];
                Real y1 = mesh.node(right_nodes[i + 1]).coords()[1];
                dy = (y1 - y0) / 2.0;
            } else if (i == right_nodes.size() - 1 && right_nodes.size() > 1) {
                Real y0 = mesh.node(right_nodes[i - 1]).coords()[1];
                Real y1 = mesh.node(right_nodes[i]).coords()[1];
                dy = (y1 - y0) / 2.0;
            }
        } else {
            // 内部节点：整个段
            Real y0 = mesh.node(right_nodes[i - 1]).coords()[1];
            Real y1 = mesh.node(right_nodes[i + 1]).coords()[1];
            dy = (y1 - y0) / 2.0;
        }
        
        Real force = shear_stress * dy * thickness;
        F_ext[node_id * dofs_per_node + 1] = force;  // Y 方向
    }
    
    LOG_INFO("Applied boundary conditions:");
    LOG_INFO("  Fixed nodes (left edge): " << right_nodes.size());
    Real total_force = 0.0;
    for (Real f : F_ext) total_force += std::abs(f);
    LOG_INFO("  Total applied force: " << total_force << " N");
    
    return F_ext;
}

int main() {
    LOG_INFO("======================================");
    LOG_INFO("Cook's Membrane - Geometric Nonlinearity");
    LOG_INFO("======================================");
    
    // ═══ 参数设置 ═══
    
    // 几何
    int nx = 8;   // X 方向单元数
    int ny = 8;   // Y 方向单元数
    
    // 材料（Neo-Hookean，近不可压）
    Real E = 250.0;         // Pa（软材料）
    Real nu = 0.4999;       // 近不可压
    
    // 载荷
    Real shear_stress = 100.0;  // Pa
    
    // 求解器参数
    NewtonRaphsonCriteria criteria;
    criteria.residual_tol = 1e-6;
    criteria.displacement_tol = 1e-8;
    criteria.max_iterations = 50;
    
    // ═══ 建立模型 ═══
    
    Model model("CooksMembrane");
    
    // 创建材料
    NeoHookean mat(E, nu, 2, true);  // 2D 平面应变
    
    // 创建网格
    int mat_id = model.add_material("NeoHookean");
    int mesh_id = model.add_mesh("membrane", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    generate_cooks_membrane_mesh(nx, ny, mesh);
    
    // ═══ 物理模块 ═══
    
    ElasticityNonlinear physics(&mat, 2, true);  // 启用几何刚度
    
    // ═══ 边界条件与载荷 ═══
    
    Assembler assembler(mesh.num_nodes(), 2);
    Vector F_ext = apply_boundary_conditions(mesh, assembler, shear_stress);
    
    // ═══ Newton-Raphson 求解 ═══
    
    LOG_INFO("\n=== Newton-Raphson Solver ===");
    
    NewtonRaphsonSolver nr_solver(model, &physics, 2);
    nr_solver.set_criteria(criteria);
    nr_solver.set_verbose(true);
    
    bool converged = nr_solver.solve(F_ext);
    
    if (!converged) {
        LOG_ERROR("Solution did not converge!");
        return 1;
    }
    
    // ═══ 结果分析 ═══
    
    const Vector& u = nr_solver.solution();
    
    // 找到右上角节点
    Index top_right_node = 0;
    Real max_x = -1.0;
    Real max_y = -1.0;
    
    for (Index node_id = 0; node_id < mesh.num_nodes(); node_id++) {
        const auto& coords = mesh.node(node_id).coords();
        
        if (coords[0] >= max_x && coords[1] >= max_y) {
            max_x = coords[0];
            max_y = coords[1];
            top_right_node = node_id;
        }
    }
    
    Real u_tip = u[top_right_node * 2 + 0];  // X 位移
    Real v_tip = u[top_right_node * 2 + 1];  // Y 位移
    
    LOG_INFO("\n=== Results ===");
    LOG_INFO("Top-right corner displacement:");
    LOG_INFO("  X displacement: " << std::fixed << std::setprecision(4) << u_tip << " mm");
    LOG_INFO("  Y displacement: " << std::fixed << std::setprecision(4) << v_tip << " mm");
    LOG_INFO("\nReference solution (Wriggers 2008):");
    LOG_INFO("  Y displacement: ~24.05 mm (geometric nonlinear)");
    LOG_INFO("  Y displacement: ~23.96 mm (linear)");
    
    Real error_percent = std::abs(v_tip - 24.05) / 24.05 * 100.0;
    LOG_INFO("\nError vs. reference: " << std::setprecision(2) << error_percent << "%");
    
    // 迭代历史
    LOG_INFO("\n=== Convergence History ===");
    const auto& history = nr_solver.iteration_history();
    for (const auto& iter : history) {
        LOG_INFO("Iteration " << iter.iteration 
                 << ": |R| = " << std::scientific << iter.residual_norm
                 << ", |du| = " << iter.displacement_norm);
    }
    
    // ═══ 输出 VTK ═══
    
    io::VTKWriter writer;
    writer.write_mesh("cooks_membrane_mesh.vtk", mesh);
    
    // 创建变形后的网格
    Mesh deformed_mesh = mesh;
    for (Index node_id = 0; node_id < mesh.num_nodes(); node_id++) {
        Vec3 coords = mesh.node(node_id).coords();
        coords[0] += u[node_id * 2 + 0];
        coords[1] += u[node_id * 2 + 1];
        deformed_mesh.node(node_id).set_coords(coords);
    }
    
    writer.write_mesh("cooks_membrane_deformed.vtk", deformed_mesh);
    
    LOG_INFO("\n=== Output ===");
    LOG_INFO("Original mesh: cooks_membrane_mesh.vtk");
    LOG_INFO("Deformed mesh: cooks_membrane_deformed.vtk");
    
    return 0;
}
