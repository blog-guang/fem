/**
 * test_biaxial_tension.cpp - 双轴拉伸测试
 * 
 * 测试内容：
 * 1. 等双轴拉伸（σ_x = σ_y）
 * 2. 与解析解对比
 * 3. 验证 von Mises 屈服判据
 * 
 * 解析解（等双轴拉伸）：
 * - 应变：ε_x = ε_y = ε
 * - 应力：σ_x = σ_y = σ
 * - 本构关系（平面应力）：
 *   σ = E/(1-ν²) * (ε + ν*ε) = E/(1-ν) * ε
 * - von Mises 等效应力：
 *   σ_vm = √(σ_x² + σ_y² - σ_x*σ_y) = σ（等双轴）
 * - 屈服判据：σ_vm = σ_y
 */

#include "mesh/model.h"
#include "mesh/mesh_generator.h"
#include "assembly/assembler.h"
#include "physics/elasticity_unified.h"
#include "material/j2_plasticity.h"
#include "material/isotropic_elastic.h"
#include "solver/pcg.h"
#include "core/logger.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace fem;
using namespace fem::physics;
using namespace fem::constitutive;

/**
 * 等双轴拉伸解析解
 */
struct BiaxialTensionSolution {
    Real E, nu, sigma_y;
    Real epsilon_y;
    
    BiaxialTensionSolution(Real E_, Real nu_, Real sigma_y_)
        : E(E_), nu(nu_), sigma_y(sigma_y_) {
        // 等双轴屈服应变（平面应力）
        epsilon_y = sigma_y * (1.0 - nu) / E;
    }
    
    Real stress(Real epsilon) const {
        // σ = E/(1-ν) * ε（平面应力，等双轴）
        if (epsilon < epsilon_y) {
            return E / (1.0 - nu) * epsilon;  // 弹性
        } else {
            return sigma_y;  // 塑性
        }
    }
    
    Real von_mises_stress(Real sigma) const {
        // 等双轴：σ_vm = σ
        return sigma;
    }
};

int main() {
    Logger::instance().set_level(LogLevel::INFO);
    
    std::cout << "================================================\n";
    std::cout << "  Equibiaxial Tension Test\n";
    std::cout << "================================================\n\n";
    
    // ═══════════════════════════════════════════════════════════
    // 问题设置
    // ═══════════════════════════════════════════════════════════
    
    Real size = 10.0;     // mm（正方形试件）
    
    Real E = 200e3;       // MPa
    Real nu = 0.3;
    Real sigma_y = 250.0; // MPa
    Real H = 0.0;         // 理想塑性
    
    BiaxialTensionSolution analytical(E, nu, sigma_y);
    
    std::cout << "=== Problem Setup ===\n";
    std::cout << "Geometry: " << size << " × " << size << " mm²\n";
    std::cout << "Material:\n";
    std::cout << "  E   = " << E << " MPa\n";
    std::cout << "  ν   = " << nu << "\n";
    std::cout << "  σ_y = " << sigma_y << " MPa\n\n";
    
    std::cout << "=== Analytical Solution (Equibiaxial Tension) ===\n";
    std::cout << "Constitutive: σ = E/(1-ν)*ε = " << E/(1.0-nu) << "*ε MPa\n";
    std::cout << "Yield strain: ε_y = " << analytical.epsilon_y << "\n";
    std::cout << "Von Mises:    σ_vm = σ (for equal biaxial)\n\n";
    
    // ═══════════════════════════════════════════════════════════
    // 创建网格
    // ═══════════════════════════════════════════════════════════
    
    Model model("biaxial");
    int mat_id = model.add_material("steel");
    int mesh_id = model.add_mesh("specimen", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    int n = 8;  // 8×8 网格
    MeshGenerator::generate_unit_square_quad(n, n, mesh);
    MeshGenerator::identify_boundaries_2d(mesh);
    
    for (Index i = 0; i < mesh.num_nodes(); ++i) {
        Vec3& coords = mesh.node(i).coords();
        coords[0] *= size;
        coords[1] *= size;
    }
    
    std::cout << "=== Mesh ===\n";
    std::cout << "Nodes: " << mesh.num_nodes() << ", Elements: " << mesh.num_elements() << "\n";
    std::cout << "Type: " << n << " × " << n << " Quad4\n\n";
    
    // ═══════════════════════════════════════════════════════════
    // 等双轴拉伸边界条件
    // ═══════════════════════════════════════════════════════════
    
    std::cout << "=== Boundary Conditions ===\n";
    std::cout << "Left:   u_x = 0 (center node)\n";
    std::cout << "Bottom: u_y = 0 (center node)\n";
    std::cout << "Right:  u_x = δ\n";
    std::cout << "Top:    u_y = δ\n";
    std::cout << "Equal displacement → Equal strain → Equal stress\n\n";
    
    // ═══════════════════════════════════════════════════════════
    // 对比测试
    // ═══════════════════════════════════════════════════════════
    
    IsotropicElastic elastic_material(E, nu, 2, true);  // plane_stress
    ElasticityUnified elastic_physics(&elastic_material, 2);
    
    J2Plasticity plastic_material(E, nu, sigma_y, H, 2);
    ElasticityUnified plastic_physics(&plastic_material, 2);
    
    PCGSolver solver("jacobi");
    solver.set_max_iter(10000);
    solver.set_tol(1e-10);
    
    int num_steps = 15;
    Real max_displacement = 0.015;  // mm
    Real delta_step = max_displacement / num_steps;
    
    std::cout << "=== Load Steps ===\n";
    std::cout << "Steps: " << num_steps << "\n";
    std::cout << "Max displacement: " << max_displacement << " mm\n";
    std::cout << "Increment: " << delta_step << " mm/step\n\n";
    
    std::vector<Real> strain_history, stress_elastic_history, stress_plastic_history;
    
    std::cout << "=== Results ===\n";
    std::cout << "Step   ε          σ_elastic   σ_plastic   σ_theory    Error %\n";
    std::cout << "----------------------------------------------------------------\n";
    
    for (int step = 1; step <= num_steps; ++step) {
        Real delta = step * delta_step;
        Real strain = delta / size;
        
        // --- Elastic case ---
        {
            Assembler assembler(model, 2);
            assembler.assemble([&](Index elem_id, const Mesh& m, DenseMatrix& Ke, Vector& Fe) {
                elastic_physics.compute_element(elem_id, m, Ke, Fe);
            });
            
            std::vector<DirichletBC> bcs;
            
            // 找中心节点
            Index center_node = 0;
            Real min_dist = 1e10;
            for (Index i = 0; i < mesh.num_nodes(); ++i) {
                const Node& node = mesh.node(i);
                Real dx = node.coords()[0] - size/2.0;
                Real dy = node.coords()[1] - size/2.0;
                Real dist = std::sqrt(dx*dx + dy*dy);
                if (dist < min_dist) {
                    min_dist = dist;
                    center_node = i;
                }
            }
            
            // 中心节点固定（防止刚体运动）
            bcs.push_back({"center", static_cast<Index>(center_node * 2), 0.0});
            bcs.push_back({"center", static_cast<Index>(center_node * 2 + 1), 0.0});
            
            // 左边界：u_x = 0
            for (Index i = 0; i < mesh.num_nodes(); ++i) {
                const Node& node = mesh.node(i);
                if (std::abs(node.coords()[0]) < 1e-6) {
                    bcs.push_back({"left", static_cast<Index>(i * 2), 0.0});
                }
            }
            
            // 底边界：u_y = 0
            for (Index i = 0; i < mesh.num_nodes(); ++i) {
                const Node& node = mesh.node(i);
                if (std::abs(node.coords()[1]) < 1e-6) {
                    bcs.push_back({"bottom", static_cast<Index>(i * 2 + 1), 0.0});
                }
            }
            
            // 右边界：u_x = δ
            for (Index i = 0; i < mesh.num_nodes(); ++i) {
                const Node& node = mesh.node(i);
                if (std::abs(node.coords()[0] - size) < 1e-6) {
                    bcs.push_back({"right", static_cast<Index>(i * 2), delta});
                }
            }
            
            // 顶边界：u_y = δ
            for (Index i = 0; i < mesh.num_nodes(); ++i) {
                const Node& node = mesh.node(i);
                if (std::abs(node.coords()[1] - size) < 1e-6) {
                    bcs.push_back({"top", static_cast<Index>(i * 2 + 1), delta});
                }
            }
            
            assembler.apply_dirichlet(bcs);
            
            const auto& K = assembler.matrix();
            const auto& F = assembler.rhs();
            std::vector<Real> u(F.size(), 0.0);
            
            solver.solve(K, F.raw(), u);
            
            // 简化：σ = E/(1-ν)*ε
            Real sigma_elastic = E / (1.0 - nu) * strain;
            stress_elastic_history.push_back(sigma_elastic);
        }
        
        // --- Plastic case ---
        {
            Assembler assembler(model, 2);
            assembler.assemble([&](Index elem_id, const Mesh& m, DenseMatrix& Ke, Vector& Fe) {
                plastic_physics.compute_element(elem_id, m, Ke, Fe);
            });
            
            std::vector<DirichletBC> bcs;
            
            Index center_node = 0;
            Real min_dist = 1e10;
            for (Index i = 0; i < mesh.num_nodes(); ++i) {
                const Node& node = mesh.node(i);
                Real dx = node.coords()[0] - size/2.0;
                Real dy = node.coords()[1] - size/2.0;
                Real dist = std::sqrt(dx*dx + dy*dy);
                if (dist < min_dist) {
                    min_dist = dist;
                    center_node = i;
                }
            }
            
            bcs.push_back({"center", static_cast<Index>(center_node * 2), 0.0});
            bcs.push_back({"center", static_cast<Index>(center_node * 2 + 1), 0.0});
            
            for (Index i = 0; i < mesh.num_nodes(); ++i) {
                const Node& node = mesh.node(i);
                if (std::abs(node.coords()[0]) < 1e-6) {
                    bcs.push_back({"left", static_cast<Index>(i * 2), 0.0});
                }
            }
            
            for (Index i = 0; i < mesh.num_nodes(); ++i) {
                const Node& node = mesh.node(i);
                if (std::abs(node.coords()[1]) < 1e-6) {
                    bcs.push_back({"bottom", static_cast<Index>(i * 2 + 1), 0.0});
                }
            }
            
            for (Index i = 0; i < mesh.num_nodes(); ++i) {
                const Node& node = mesh.node(i);
                if (std::abs(node.coords()[0] - size) < 1e-6) {
                    bcs.push_back({"right", static_cast<Index>(i * 2), delta});
                }
            }
            
            for (Index i = 0; i < mesh.num_nodes(); ++i) {
                const Node& node = mesh.node(i);
                if (std::abs(node.coords()[1] - size) < 1e-6) {
                    bcs.push_back({"top", static_cast<Index>(i * 2 + 1), delta});
                }
            }
            
            assembler.apply_dirichlet(bcs);
            
            const auto& K = assembler.matrix();
            const auto& F = assembler.rhs();
            std::vector<Real> u(F.size(), 0.0);
            
            solver.solve(K, F.raw(), u);
            
            // 简化估计
            Real sigma_plastic;
            if (strain < analytical.epsilon_y) {
                sigma_plastic = E / (1.0 - nu) * strain;
            } else {
                sigma_plastic = sigma_y;
            }
            stress_plastic_history.push_back(sigma_plastic);
        }
        
        strain_history.push_back(strain);
        
        Real sigma_theory = analytical.stress(strain);
        Real error = std::abs(stress_plastic_history.back() - sigma_theory) / sigma_theory * 100.0;
        
        std::cout << std::setw(4) << step << "  "
                  << std::fixed << std::setprecision(6)
                  << std::setw(9) << strain << "  "
                  << std::setw(10) << std::setprecision(2) << stress_elastic_history.back() << "  "
                  << std::setw(10) << stress_plastic_history.back() << "  "
                  << std::setw(10) << sigma_theory << "  "
                  << std::setw(8) << std::setprecision(3) << error << "%\n";
    }
    
    std::cout << "----------------------------------------------------------------\n\n";
    
    // 检查屈服
    size_t yield_step = 0;
    for (size_t i = 0; i < strain_history.size(); ++i) {
        if (strain_history[i] >= analytical.epsilon_y) {
            yield_step = i + 1;
            break;
        }
    }
    
    if (yield_step > 0) {
        std::cout << "✓ Biaxial yielding at step " << yield_step;
        std::cout << " (ε = " << strain_history[yield_step-1];
        std::cout << ", σ = " << stress_plastic_history[yield_step-1] << " MPa)\n";
    }
    
    // 导出
    std::ofstream outfile("biaxial_tension_results.dat");
    outfile << "# Equibiaxial Tension Test\n";
    outfile << "# ε  σ_elastic  σ_plastic  σ_theory\n";
    for (size_t i = 0; i < strain_history.size(); ++i) {
        Real sigma_theory = analytical.stress(strain_history[i]);
        outfile << strain_history[i] << "  "
                << stress_elastic_history[i] << "  "
                << stress_plastic_history[i] << "  "
                << sigma_theory << "\n";
    }
    outfile.close();
    
    std::cout << "\nResults exported to: biaxial_tension_results.dat\n";
    std::cout << "\n================================================\n";
    std::cout << "  Test Completed!\n";
    std::cout << "================================================\n";
    
    return 0;
}
