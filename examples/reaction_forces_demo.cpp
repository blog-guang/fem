/**
 * test_reaction_forces.cpp - 反力提取验证
 * 
 * 测试场景：
 * 1. 单轴拉伸杆件
 * 2. 左端固定 (u=0)
 * 3. 右端施加位移 (u=δ)
 * 4. 验证反力 R = -F_applied
 */

#include "mesh/model.h"
#include "mesh/mesh_generator.h"
#include "assembly/assembler.h"
#include "physics/elasticity_unified.h"
#include "material/isotropic_elastic.h"
#include "math/cg.h"
#include "core/logger.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace fem;
using namespace fem::physics;
using namespace fem::constitutive;

int main() {
    std::cout << "============================================\n";
    std::cout << "  Reaction Force Extraction Test\n";
    std::cout << "============================================\n\n";

    Logger::instance().set_level(LogLevel::WARN);

    // ========== 参数设置 ==========
    const Real length = 10.0;      // 长度 (mm)
    const Real width = 1.0;        // 宽度 (mm)
    const Real E = 200000.0;       // 杨氏模量 (MPa)
    const Real nu = 0.3;           // 泊松比
    const Real delta = 0.01;       // 施加位移 (mm)
    
    std::cout << "=== Problem Setup ===\n";
    std::cout << "Geometry: " << length << " × " << width << " mm²\n";
    std::cout << "Material: E = " << E << " MPa, ν = " << nu << "\n";
    std::cout << "Applied displacement: δ = " << delta << " mm\n\n";

    // ========== 解析解 ==========
    // 平面应力状态，考虑泊松效应
    const Real strain = delta / length;
    const Real stress = E / (1.0 - nu * nu) * strain;  // σ = E/(1-ν²) · ε
    const Real force_analytical = stress * width;  // F = σ * A
    
    std::cout << "=== Analytical Solution (Plane Stress) ===\n";
    std::cout << "Strain: ε_x = δ/L = " << strain << "\n";
    std::cout << "Stress: σ_x = E/(1-ν²)·ε = " << stress << " MPa\n";
    std::cout << "Reaction force: F = σ*A = " << force_analytical << " N\n\n";

    // ========== 创建网格 ==========
    Model model("reaction_test");
    int mat_id = model.add_material("steel");
    int mesh_id = model.add_mesh("bar", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    MeshGenerator::generate_unit_square_quad(5, 1, mesh);  // 5x1 网格
    MeshGenerator::identify_boundaries_2d(mesh);
    
    // 缩放到实际长度
    for (Index i = 0; i < mesh.num_nodes(); ++i) {
        Vec3& coords = mesh.node(i).coords();
        coords[0] *= length;  // 只缩放 x 方向
        coords[1] *= width;   // y 方向保持 width
    }
    
    std::cout << "=== Mesh ===\n";
    std::cout << "Nodes: " << mesh.num_nodes() << ", Elements: " << mesh.num_elements() << "\n";
    std::cout << "Type: 5 × 1 Quad4\n\n";

    // ========== 材料和物理场 ==========
    IsotropicElastic material(E, nu, 2, true);  // 2D plane stress
    ElasticityUnified physics(&material, 2);

    // ========== 装配系统 ==========
    Assembler assembler(model, 2);
    assembler.assemble([&](Index elem_id, const Mesh& m, DenseMatrix& Ke, Vector& Fe) {
        physics.compute_element(elem_id, m, Ke, Fe);
    });

    // ========== 边界条件 ==========
    std::vector<DirichletBC> bcs;
    
    // 左端固定: u_x = 0, u_y = 0
    bcs.push_back({"left", 0, 0.0});
    bcs.push_back({"left", 1, 0.0});
    
    // 右端施加位移: u_x = δ
    bcs.push_back({"right", 0, delta});
    
    // 底边和顶边 u_y = 0 (防止刚体位移)
    bcs.push_back({"bottom", 1, 0.0});
    bcs.push_back({"top", 1, 0.0});

    assembler.apply_dirichlet(bcs);

    std::cout << "=== Boundary Conditions ===\n";
    std::cout << "Left:   u_x = u_y = 0 (fixed)\n";
    std::cout << "Right:  u_x = " << delta << " mm (applied)\n";
    std::cout << "Top/Bottom: u_y = 0 (constraint)\n\n";

    // ========== 求解 ==========
    const auto& K = assembler.matrix();
    const auto& F = assembler.rhs();
    Vector u(F.size(), 0.0);

    CGSolver cg_solver;
    cg_solver.set_tol(1e-8);
    cg_solver.set_max_iter(1000);
    auto result = cg_solver.solve(K, F, u);

    if (!result.converged) {
        std::cerr << "❌ CG solver failed!\n";
        return 1;
    }

    std::cout << "=== Solution ===\n";
    std::cout << "CG iterations: " << result.iterations << "\n";
    std::cout << "Residual: " << result.residual << "\n\n";

    // ========== 提取反力 ==========
    Vector R = assembler.compute_reaction_forces(u.raw());

    // 计算左端节点的总反力 (x 方向)
    Real total_reaction_x = 0.0;
    Real total_reaction_y = 0.0;
    
    const auto& left_nodes = mesh.boundary("left");
    for (Index node_id : left_nodes) {
        Index dof_x = node_id * 2;
        Index dof_y = node_id * 2 + 1;
        total_reaction_x += R[dof_x];
        total_reaction_y += R[dof_y];
    }

    std::cout << "=== Reaction Forces (Left Boundary) ===\n";
    std::cout << "Total R_x: " << total_reaction_x << " N\n";
    std::cout << "Total R_y: " << total_reaction_y << " N\n\n";

    // ========== 验证 ==========
    Real error_x = std::abs(total_reaction_x + force_analytical) / std::abs(force_analytical) * 100.0;
    Real error_y = std::abs(total_reaction_y);

    std::cout << "=== Verification ===\n";
    std::cout << "Expected: R_x = -" << force_analytical << " N (opposite to applied force)\n";
    std::cout << "FEM:      R_x = " << total_reaction_x << " N\n";
    std::cout << "Error:    " << std::fixed << std::setprecision(3) << error_x << " %\n\n";

    std::cout << "Expected: R_y ≈ 0 N\n";
    std::cout << "FEM:      R_y = " << std::scientific << std::setprecision(3) << total_reaction_y << " N\n\n";

    // ========== 详细输出 ==========
    std::cout << "=== Detailed Reaction Forces (Left Nodes) ===\n";
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Node    R_x (N)     R_y (N)\n";
    std::cout << "--------------------------------\n";
    
    for (Index node_id : left_nodes) {
        Index dof_x = node_id * 2;
        Index dof_y = node_id * 2 + 1;
        std::cout << std::setw(4) << node_id << "  "
                  << std::setw(10) << R[dof_x] << "  "
                  << std::setw(10) << R[dof_y] << "\n";
    }
    std::cout << "--------------------------------\n";
    std::cout << "Sum   " << std::setw(10) << total_reaction_x << "  "
              << std::setw(10) << total_reaction_y << "\n\n";

    // ========== 结果判断 ==========
    if (error_x < 1.0 && error_y < 1e-6) {
        std::cout << "✅ Test PASSED!\n";
        std::cout << "   Reaction forces match analytical solution.\n";
    } else {
        std::cout << "❌ Test FAILED!\n";
        std::cout << "   Error too large.\n";
    }

    std::cout << "\n============================================\n";
    std::cout << "  Test Completed\n";
    std::cout << "============================================\n";

    return 0;
}
