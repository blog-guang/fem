/**
 * test_plasticity_simple.cpp - 简单拉伸塑性测试
 * 
 * 问题描述：
 * - 单轴拉伸试件
 * - J2 塑性材料（理想塑性，H=0）
 * - 与解析解对比验证
 * 
 * 解析解（单轴拉伸）：
 * - 弹性阶段：σ = E*ε
 * - 塑性阶段：σ = σ_y (理想塑性)
 * - 应变：ε = ε_e + ε_p
 * 
 * 这是最简单的塑性问题，有精确解析解。
 */

#include "mesh/model.h"
#include "mesh/mesh_generator.h"
#include "assembly/assembler.h"
#include "physics/elasticity_unified.h"
#include "material/j2_plasticity.h"
#include "material/isotropic_elastic.h"
#include "solver/cg.h"
#include "core/logger.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace fem;
using namespace fem::physics;
using namespace fem::constitutive;

int main() {
    Logger::instance().set_level(LogLevel::INFO);
    
    std::cout << "========================================\n";
    std::cout << "  Simple Uniaxial Tension Plasticity\n";
    std::cout << "========================================\n\n";
    
    // ═══════════════════════════════════════════════════════════
    // 几何和材料参数
    // ═══════════════════════════════════════════════════════════
    
    // 试件几何
    Real length = 10.0;   // 长度
    Real width = 1.0;     // 宽度
    
    // 材料参数
    Real E = 200e3;        // MPa（钢）
    Real nu = 0.3;
    Real sigma_y = 250.0;  // MPa（屈服应力）
    Real H = 0.0;          // 理想塑性（无硬化）
    
    std::cout << "=== Problem Setup ===\n";
    std::cout << "Geometry:\n";
    std::cout << "  Length L = " << length << " m\n";
    std::cout << "  Width  b = " << width << " m\n";
    std::cout << "\nMaterial (J2 Plasticity, Perfect Plastic):\n";
    std::cout << "  E   = " << E << " MPa\n";
    std::cout << "  ν   = " << nu << "\n";
    std::cout << "  σ_y = " << sigma_y << " MPa\n";
    std::cout << "  H   = " << H << " MPa (perfect plastic)\n\n";
    
    // 解析解
    Real epsilon_y = sigma_y / E;  // 屈服应变
    
    std::cout << "=== Analytical Solution ===\n";
    std::cout << "Yield strain: ε_y = σ_y/E = " << epsilon_y << "\n";
    std::cout << "Elastic phase:  σ = E*ε     (ε < ε_y)\n";
    std::cout << "Plastic phase:  σ = σ_y     (ε > ε_y)\n\n";
    
    // ═══════════════════════════════════════════════════════════
    // 创建网格
    // ═══════════════════════════════════════════════════════════
    
    Model model("uniaxial_tension");
    int mat_id = model.add_material("steel");
    int mesh_id = model.add_mesh("specimen", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    // 生成矩形网格
    int nx = 10;
    int ny = 2;
    
    MeshGenerator::generate_unit_square_quad(nx, ny, mesh);
    
    // 缩放到实际尺寸
    for (Index i = 0; i < mesh.num_nodes(); ++i) {
        Node& node = mesh.node(i);
        Vec3& coords = node.coords();
        coords[0] = coords[0] * length;  // x: [0, L]
        coords[1] = coords[1] * width;   // y: [0, b]
    }
    
    std::cout << "=== Mesh ===\n";
    std::cout << "  Nodes:    " << mesh.num_nodes() << "\n";
    std::cout << "  Elements: " << mesh.num_elements() << "\n";
    std::cout << "  Mesh:     " << nx << " × " << ny << " Quad4\n\n";
    
    // ═══════════════════════════════════════════════════════════
    // 对比测试：弹性 vs 塑性
    // ═══════════════════════════════════════════════════════════
    
    std::cout << "=== Running Two Cases ===\n\n";
    
    // Case 1: 弹性材料
    std::cout << "Case 1: Elastic Material\n";
    std::cout << "------------------------\n";
    
    IsotropicElastic elastic_material(E, nu, 2, true);  // 2D, plane_stress
    ElasticityUnified elastic_physics(&elastic_material, 2);
    
    std::vector<Real> strain_elastic, stress_elastic;
    
    // Case 2: 塑性材料
    std::cout << "\nCase 2: Plastic Material (J2, H=0)\n";
    std::cout << "-----------------------------------\n";
    
    J2Plasticity plastic_material(E, nu, sigma_y, H, 2);  // 2D
    ElasticityUnified plastic_physics(&plastic_material, 2);
    
    std::vector<Real> strain_plastic, stress_plastic;
    
    // ═══════════════════════════════════════════════════════════
    // 增量加载
    // ═══════════════════════════════════════════════════════════
    
    int num_steps = 20;
    Real max_displacement = 0.02;  // 最大位移（2% 应变）
    Real du = max_displacement / num_steps;
    
    std::cout << "\nLoad Steps: " << num_steps << "\n";
    std::cout << "Max displacement: " << max_displacement << " m (" 
              << (max_displacement/length*100) << "% strain)\n\n";
    
    std::cout << "Step  Displacement   Elastic σ   Plastic σ   σ_plastic/σ_y\n";
    std::cout << "-----------------------------------------------------------\n";
    
    for (int step = 1; step <= num_steps; ++step) {
        Real u_applied = step * du;
        Real strain_applied = u_applied / length;
        
        // --- Case 1: Elastic ---
        {
            Assembler assembler(model, 2);
            
            assembler.assemble([&](Index elem_id, const Mesh& m, DenseMatrix& Ke, Vector& Fe) {
                elastic_physics.compute_element(elem_id, m, Ke, Fe);
            });
            
            // 边界条件
            std::vector<DirichletBC> bcs;
            
            // 左端固定
            for (Index i = 0; i < mesh.num_nodes(); ++i) {
                const Node& node = mesh.node(i);
                if (std::abs(node.coords()[0]) < 1e-6) {
                    bcs.push_back({"left", static_cast<Index>(i * 2), 0.0});  // u_x = 0
                }
            }
            
            // 右端施加位移
            for (Index i = 0; i < mesh.num_nodes(); ++i) {
                const Node& node = mesh.node(i);
                if (std::abs(node.coords()[0] - length) < 1e-6) {
                    bcs.push_back({"right", static_cast<Index>(i * 2), u_applied});  // u_x = u
                }
            }
            
            // 底部固定 y 方向（防止刚体运动）
            for (Index i = 0; i < mesh.num_nodes(); ++i) {
                const Node& node = mesh.node(i);
                if (std::abs(node.coords()[1]) < 1e-6) {
                    bcs.push_back({"bottom", static_cast<Index>(i * 2 + 1), 0.0});  // u_y = 0
                }
            }
            
            assembler.apply_dirichlet(bcs);
            
            // 求解
            const auto& K = assembler.matrix();
            const auto& F = assembler.rhs();
            std::vector<Real> u(F.size(), 0.0);
            
            CGSolver cg;
            cg.set_tol(1e-8);
            cg.solve(K, F.raw(), u);
            
            // 提取应力（简化：从反力计算）
            // σ = F_reaction / A
            Real F_total = 0.0;
            for (Index i = 0; i < mesh.num_nodes(); ++i) {
                const Node& node = mesh.node(i);
                if (std::abs(node.coords()[0]) < 1e-6) {
                    F_total += std::abs((K * u)[i * 2]);  // 反力
                }
            }
            
            Real area = width * 1.0;  // 单位厚度
            Real sigma_elastic = F_total / area;
            
            strain_elastic.push_back(strain_applied);
            stress_elastic.push_back(sigma_elastic);
        }
        
        // --- Case 2: Plastic ---
        {
            Assembler assembler(model, 2);
            
            assembler.assemble([&](Index elem_id, const Mesh& m, DenseMatrix& Ke, Vector& Fe) {
                plastic_physics.compute_element(elem_id, m, Ke, Fe);
            });
            
            // 边界条件（同上）
            std::vector<DirichletBC> bcs;
            
            for (Index i = 0; i < mesh.num_nodes(); ++i) {
                const Node& node = mesh.node(i);
                if (std::abs(node.coords()[0]) < 1e-6) {
                    bcs.push_back({"left", static_cast<Index>(i * 2), 0.0});
                }
            }
            
            for (Index i = 0; i < mesh.num_nodes(); ++i) {
                const Node& node = mesh.node(i);
                if (std::abs(node.coords()[0] - length) < 1e-6) {
                    bcs.push_back({"right", static_cast<Index>(i * 2), u_applied});
                }
            }
            
            for (Index i = 0; i < mesh.num_nodes(); ++i) {
                const Node& node = mesh.node(i);
                if (std::abs(node.coords()[1]) < 1e-6) {
                    bcs.push_back({"bottom", static_cast<Index>(i * 2 + 1), 0.0});
                }
            }
            
            assembler.apply_dirichlet(bcs);
            
            // 求解
            const auto& K = assembler.matrix();
            const auto& F = assembler.rhs();
            std::vector<Real> u(F.size(), 0.0);
            
            CGSolver cg;
            cg.set_tol(1e-8);
            cg.solve(K, F.raw(), u);
            
            // 提取应力
            Real F_total = 0.0;
            for (Index i = 0; i < mesh.num_nodes(); ++i) {
                const Node& node = mesh.node(i);
                if (std::abs(node.coords()[0]) < 1e-6) {
                    F_total += std::abs((K * u)[i * 2]);
                }
            }
            
            Real area = width * 1.0;
            Real sigma_plastic = F_total / area;
            
            strain_plastic.push_back(strain_applied);
            stress_plastic.push_back(sigma_plastic);
        }
        
        // 打印结果
        printf("%3d   %10.6f    %9.2f    %9.2f      %6.3f\n",
               step,
               u_applied,
               stress_elastic.back(),
               stress_plastic.back(),
               stress_plastic.back() / sigma_y);
    }
    
    std::cout << "-----------------------------------------------------------\n\n";
    
    // ═══════════════════════════════════════════════════════════
    // 结果对比
    // ═══════════════════════════════════════════════════════════
    
    std::cout << "=== Results Summary ===\n\n";
    
    // 找到屈服点
    size_t yield_index = 0;
    for (size_t i = 0; i < stress_plastic.size(); ++i) {
        if (stress_plastic[i] >= 0.99 * sigma_y) {
            yield_index = i;
            break;
        }
    }
    
    if (yield_index > 0) {
        std::cout << "Yielding occurred at:\n";
        std::cout << "  Strain:  ε = " << strain_plastic[yield_index] << "\n";
        std::cout << "  Stress:  σ = " << stress_plastic[yield_index] << " MPa\n";
        std::cout << "  σ/σ_y = " << stress_plastic[yield_index]/sigma_y << "\n\n";
        
        std::cout << "Analytical yield strain: ε_y = " << epsilon_y << "\n";
        std::cout << "FEM yield strain:        ε   = " << strain_plastic[yield_index] << "\n";
        std::cout << "Relative error: " 
                  << std::abs(strain_plastic[yield_index] - epsilon_y)/epsilon_y * 100 
                  << "%\n\n";
    }
    
    // ═══════════════════════════════════════════════════════════
    // 导出结果
    // ═══════════════════════════════════════════════════════════
    
    std::ofstream outfile("plasticity_simple_results.dat");
    outfile << "# Simple Uniaxial Tension Plasticity\n";
    outfile << "# ε           σ_elastic    σ_plastic    σ_analytical\n";
    
    for (size_t i = 0; i < strain_elastic.size(); ++i) {
        Real sigma_analytical = (strain_elastic[i] < epsilon_y) ? 
                                E * strain_elastic[i] : sigma_y;
        
        outfile << strain_elastic[i] << "  "
                << stress_elastic[i] << "  "
                << stress_plastic[i] << "  "
                << sigma_analytical << "\n";
    }
    
    outfile.close();
    std::cout << "Results exported to: plasticity_simple_results.dat\n\n";
    
    std::cout << "========================================\n";
    std::cout << "  Test Completed Successfully!\n";
    std::cout << "========================================\n";
    
    return 0;
}
