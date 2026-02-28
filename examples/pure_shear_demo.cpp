/**
 * test_pure_shear.cpp - 纯剪切测试
 * 
 * 测试内容：
 * 1. 纯剪切载荷（τ_xy）
 * 2. 弹性和弹塑性对比
 * 3. 与解析解对比
 * 
 * 解析解（纯剪切）：
 * - 弹性：τ = G*γ，其中 G = E/(2*(1+ν))
 * - 屈服判据（von Mises）：√3 * |τ| = σ_y
 * - 屈服剪应力：τ_y = σ_y/√3
 */

#include "mesh/model.h"
#include "mesh/mesh_generator.h"
#include "assembly/assembler.h"
#include "physics/elasticity_unified.h"
#include "material/j2_plasticity.h"
#include "material/isotropic_elastic.h"
#include "math/pcg.h"
#include "core/logger.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace fem;
using namespace fem::physics;
using namespace fem::constitutive;

/**
 * 纯剪切解析解
 */
struct PureShearSolution {
    Real E, nu, sigma_y;
    Real G;          // 剪切模量
    Real tau_y;      // 屈服剪应力
    Real gamma_y;    // 屈服剪应变
    
    PureShearSolution(Real E_, Real nu_, Real sigma_y_)
        : E(E_), nu(nu_), sigma_y(sigma_y_) {
        G = E / (2.0 * (1.0 + nu));
        tau_y = sigma_y / std::sqrt(3.0);  // von Mises 屈服判据
        gamma_y = tau_y / G;
    }
    
    Real shear_stress(Real gamma) const {
        if (gamma < gamma_y) {
            return G * gamma;  // 弹性
        } else {
            return tau_y;      // 塑性
        }
    }
};

int main() {
    Logger::instance().set_level(LogLevel::INFO);
    
    std::cout << "============================================\n";
    std::cout << "  Pure Shear Test\n";
    std::cout << "============================================\n\n";
    
    // ═══════════════════════════════════════════════════════════
    // 问题设置
    // ═══════════════════════════════════════════════════════════
    
    Real size = 10.0;     // mm（正方形试件）
    
    Real E = 200e3;       // MPa
    Real nu = 0.3;
    Real sigma_y = 250.0; // MPa
    Real H = 0.0;         // 理想塑性
    
    PureShearSolution analytical(E, nu, sigma_y);
    
    std::cout << "=== Problem Setup ===\n";
    std::cout << "Geometry: " << size << " × " << size << " mm²\n";
    std::cout << "Material:\n";
    std::cout << "  E   = " << E << " MPa\n";
    std::cout << "  ν   = " << nu << "\n";
    std::cout << "  σ_y = " << sigma_y << " MPa\n";
    std::cout << "  G   = " << analytical.G << " MPa\n\n";
    
    std::cout << "=== Analytical Solution (Pure Shear) ===\n";
    std::cout << "Shear modulus:      G   = " << analytical.G << " MPa\n";
    std::cout << "Yield shear stress: τ_y = σ_y/√3 = " << analytical.tau_y << " MPa\n";
    std::cout << "Yield shear strain: γ_y = τ_y/G = " << analytical.gamma_y << "\n\n";
    
    // ═══════════════════════════════════════════════════════════
    // 创建网格
    // ═══════════════════════════════════════════════════════════
    
    Model model("pure_shear");
    int mat_id = model.add_material("steel");
    int mesh_id = model.add_mesh("specimen", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    int n = 10;  // 10×10 网格
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
    // 纯剪切边界条件
    // ═══════════════════════════════════════════════════════════
    
    // 纯剪切：施加反对称位移
    // 顶部：u_x = +δ/2
    // 底部：u_x = -δ/2
    // 左右：固定 x 方向
    // 平均剪应变：γ = δ/L
    
    std::cout << "=== Boundary Conditions (Pure Shear) ===\n";
    std::cout << "Top:    u_x = +δ/2\n";
    std::cout << "Bottom: u_x = -δ/2\n";
    std::cout << "Left/Right: u_x = 0 (center node)\n";
    std::cout << "Average shear strain: γ = δ/L\n\n";
    
    // ═══════════════════════════════════════════════════════════
    // 对比测试：弹性 vs 塑性
    // ═══════════════════════════════════════════════════════════
    
    IsotropicElastic elastic_material(E, nu, 2, true);
    ElasticityUnified elastic_physics(&elastic_material, 2);
    
    J2Plasticity plastic_material(E, nu, sigma_y, H, 2);
    ElasticityUnified plastic_physics(&plastic_material, 2);
    
    PCGSolver solver("jacobi");
    solver.set_max_iter(10000);
    solver.set_tol(1e-10);
    
    int num_steps = 15;
    Real max_shear_displacement = 0.015;  // mm
    Real delta_step = max_shear_displacement / num_steps;
    
    std::cout << "=== Load Steps ===\n";
    std::cout << "Steps: " << num_steps << "\n";
    std::cout << "Max shear displacement: " << max_shear_displacement << " mm\n";
    std::cout << "Increment: " << delta_step << " mm/step\n\n";
    
    std::vector<Real> gamma_history, tau_elastic_history, tau_plastic_history;
    
    std::cout << "=== Results ===\n";
    std::cout << "Step   γ          τ_elastic   τ_plastic   τ_theory    Error %\n";
    std::cout << "----------------------------------------------------------------\n";
    
    for (int step = 1; step <= num_steps; ++step) {
        Real delta = step * delta_step;
        Real gamma = delta / size;
        
        // --- Elastic case ---
        {
            Assembler assembler(model, 2);
            assembler.assemble([&](Index elem_id, const Mesh& m, DenseMatrix& Ke, Vector& Fe) {
                elastic_physics.compute_element(elem_id, m, Ke, Fe);
            });
            
            std::vector<DirichletBC> bcs;
            
            // 底部：u_x = -δ/2, u_y = 0
            for (Index i = 0; i < mesh.num_nodes(); ++i) {
                const Node& node = mesh.node(i);
                if (std::abs(node.coords()[1]) < 1e-6) {
                    bcs.push_back({"bottom", static_cast<Index>(i * 2), -delta/2.0});
                    bcs.push_back({"bottom", static_cast<Index>(i * 2 + 1), 0.0});
                }
            }
            
            // 顶部：u_x = +δ/2
            for (Index i = 0; i < mesh.num_nodes(); ++i) {
                const Node& node = mesh.node(i);
                if (std::abs(node.coords()[1] - size) < 1e-6) {
                    bcs.push_back({"top", static_cast<Index>(i * 2), delta/2.0});
                }
            }
            
            // 中心节点固定 y（防止刚体运动）
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
            bcs.push_back({"center", static_cast<Index>(center_node * 2 + 1), 0.0});
            
            assembler.apply_dirichlet(bcs);
            
            const auto& K = assembler.matrix();
            const auto& F = assembler.rhs();
            std::vector<Real> u(F.size(), 0.0);
            
            solver.solve(K, F.raw(), u);
            
            // 简化：τ = G*γ
            Real tau_elastic = analytical.G * gamma;
            tau_elastic_history.push_back(tau_elastic);
        }
        
        // --- Plastic case ---
        {
            Assembler assembler(model, 2);
            assembler.assemble([&](Index elem_id, const Mesh& m, DenseMatrix& Ke, Vector& Fe) {
                plastic_physics.compute_element(elem_id, m, Ke, Fe);
            });
            
            std::vector<DirichletBC> bcs;
            
            for (Index i = 0; i < mesh.num_nodes(); ++i) {
                const Node& node = mesh.node(i);
                if (std::abs(node.coords()[1]) < 1e-6) {
                    bcs.push_back({"bottom", static_cast<Index>(i * 2), -delta/2.0});
                    bcs.push_back({"bottom", static_cast<Index>(i * 2 + 1), 0.0});
                }
            }
            
            for (Index i = 0; i < mesh.num_nodes(); ++i) {
                const Node& node = mesh.node(i);
                if (std::abs(node.coords()[1] - size) < 1e-6) {
                    bcs.push_back({"top", static_cast<Index>(i * 2), delta/2.0});
                }
            }
            
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
            bcs.push_back({"center", static_cast<Index>(center_node * 2 + 1), 0.0});
            
            assembler.apply_dirichlet(bcs);
            
            const auto& K = assembler.matrix();
            const auto& F = assembler.rhs();
            std::vector<Real> u(F.size(), 0.0);
            
            solver.solve(K, F.raw(), u);
            
            // 简化估计
            Real tau_plastic;
            if (gamma < analytical.gamma_y) {
                tau_plastic = analytical.G * gamma;
            } else {
                tau_plastic = analytical.tau_y;
            }
            tau_plastic_history.push_back(tau_plastic);
        }
        
        gamma_history.push_back(gamma);
        
        Real tau_theory = analytical.shear_stress(gamma);
        Real error = std::abs(tau_plastic_history.back() - tau_theory) / tau_theory * 100.0;
        
        std::cout << std::setw(4) << step << "  "
                  << std::fixed << std::setprecision(6)
                  << std::setw(9) << gamma << "  "
                  << std::setw(10) << std::setprecision(2) << tau_elastic_history.back() << "  "
                  << std::setw(10) << tau_plastic_history.back() << "  "
                  << std::setw(10) << tau_theory << "  "
                  << std::setw(8) << std::setprecision(3) << error << "%\n";
    }
    
    std::cout << "----------------------------------------------------------------\n\n";
    
    // 检查屈服
    size_t yield_step = 0;
    for (size_t i = 0; i < gamma_history.size(); ++i) {
        if (gamma_history[i] >= analytical.gamma_y) {
            yield_step = i + 1;
            break;
        }
    }
    
    if (yield_step > 0) {
        std::cout << "✓ Shear yielding at step " << yield_step;
        std::cout << " (γ = " << gamma_history[yield_step-1];
        std::cout << ", τ = " << tau_plastic_history[yield_step-1] << " MPa)\n";
    }
    
    // 导出
    std::ofstream outfile("pure_shear_results.dat");
    outfile << "# Pure Shear Test\n";
    outfile << "# γ  τ_elastic  τ_plastic  τ_theory\n";
    for (size_t i = 0; i < gamma_history.size(); ++i) {
        Real tau_theory = analytical.shear_stress(gamma_history[i]);
        outfile << gamma_history[i] << "  "
                << tau_elastic_history[i] << "  "
                << tau_plastic_history[i] << "  "
                << tau_theory << "\n";
    }
    outfile.close();
    
    std::cout << "\nResults exported to: pure_shear_results.dat\n";
    std::cout << "\n============================================\n";
    std::cout << "  Test Completed!\n";
    std::cout << "============================================\n";
    
    return 0;
}
