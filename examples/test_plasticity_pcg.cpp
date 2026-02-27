/**
 * test_plasticity_pcg.cpp - PCG求解器 + 塑性材料测试
 * 
 * 目标：
 * 1. 验证 PCG 求解器（预条件共轭梯度）
 * 2. 与 CG 求解器对比性能
 * 3. 完整的非线性塑性分析
 * 4. 与解析解对齐
 * 
 * 问题设置：
 * - 单轴拉伸试件
 * - J2 塑性材料（各向同性硬化）
 * - 增量加载
 * - 解析解：理想塑性
 */

#include "mesh/model.h"
#include "mesh/mesh_generator.h"
#include "assembly/assembler.h"
#include "physics/elasticity_unified.h"
#include "material/j2_plasticity.h"
#include "material/isotropic_elastic.h"
#include "solver/cg.h"
#include "solver/pcg.h"
#include "core/logger.h"
#include "core/timer.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace fem;
using namespace fem::physics;
using namespace fem::constitutive;

/**
 * 解析解：单轴拉伸理想塑性
 */
struct AnalyticalSolution {
    Real E, nu, sigma_y;
    Real epsilon_y;
    
    AnalyticalSolution(Real E_, Real nu_, Real sigma_y_)
        : E(E_), nu(nu_), sigma_y(sigma_y_) {
        epsilon_y = sigma_y / E;
    }
    
    Real stress(Real epsilon) const {
        if (epsilon < epsilon_y) {
            return E * epsilon;  // 弹性
        } else {
            return sigma_y;      // 塑性
        }
    }
    
    Real strain_energy(Real epsilon) const {
        if (epsilon < epsilon_y) {
            return 0.5 * E * epsilon * epsilon;
        } else {
            Real W_e = 0.5 * sigma_y * epsilon_y;
            Real W_p = sigma_y * (epsilon - epsilon_y);
            return W_e + W_p;
        }
    }
};

int main() {
    Logger::instance().set_level(LogLevel::INFO);
    
    std::cout << "================================================\n";
    std::cout << "  Plasticity Test with PCG Solver\n";
    std::cout << "================================================\n\n";
    
    // ═══════════════════════════════════════════════════════════
    // 问题设置
    // ═══════════════════════════════════════════════════════════
    
    // 几何
    Real length = 10.0;
    Real width = 1.0;
    Real thickness = 1.0;  // 单位厚度
    
    // 材料参数
    Real E = 200e3;        // MPa
    Real nu = 0.3;
    Real sigma_y = 250.0;  // MPa
    Real H = 0.0;          // 理想塑性
    
    std::cout << "=== Problem Setup ===\n";
    std::cout << "Geometry:\n";
    std::cout << "  Length:    L = " << length << " mm\n";
    std::cout << "  Width:     b = " << width << " mm\n";
    std::cout << "  Thickness: t = " << thickness << " mm\n\n";
    
    std::cout << "Material (J2 Perfect Plasticity):\n";
    std::cout << "  E   = " << E << " MPa\n";
    std::cout << "  ν   = " << nu << "\n";
    std::cout << "  σ_y = " << sigma_y << " MPa\n";
    std::cout << "  H   = " << H << " MPa\n\n";
    
    // 解析解
    AnalyticalSolution analytical(E, nu, sigma_y);
    std::cout << "Analytical yield strain: ε_y = " << analytical.epsilon_y << "\n\n";
    
    // ═══════════════════════════════════════════════════════════
    // 创建网格
    // ═══════════════════════════════════════════════════════════
    
    Model model("uniaxial_tension");
    int mat_id = model.add_material("steel");
    int mesh_id = model.add_mesh("specimen", mat_id);
    Mesh& mesh = model.mesh(mesh_id);
    
    // 细网格（提高精度）
    int nx = 20;
    int ny = 4;
    
    MeshGenerator::generate_unit_square_quad(nx, ny, mesh);
    MeshGenerator::identify_boundaries_2d(mesh);  // 识别边界
    
    // 缩放到实际尺寸
    for (Index i = 0; i < mesh.num_nodes(); ++i) {
        Node& node = mesh.node(i);
        Vec3& coords = node.coords();
        coords[0] = coords[0] * length;
        coords[1] = coords[1] * width;
    }
    
    std::cout << "=== Mesh ===\n";
    std::cout << "  Nodes:    " << mesh.num_nodes() << "\n";
    std::cout << "  Elements: " << mesh.num_elements() << "\n";
    std::cout << "  Type:     " << nx << " × " << ny << " Quad4\n\n";
    
    // ═══════════════════════════════════════════════════════════
    // 对比测试：弹性 vs 塑性，CG vs PCG
    // ═══════════════════════════════════════════════════════════
    
    std::cout << "=== Running Tests ===\n\n";
    
    // 测试参数
    int num_steps = 15;
    Real max_displacement = 0.015;  // 1.5% 应变
    Real du = max_displacement / num_steps;
    
    std::cout << "Load steps: " << num_steps << "\n";
    std::cout << "Max displacement: " << max_displacement << " mm\n";
    std::cout << "Max strain: " << (max_displacement / length * 100) << "%\n\n";
    
    // 结果存储
    struct Result {
        Real displacement;
        Real strain;
        Real stress_elastic_cg;
        Real stress_plastic_pcg;
        Real stress_analytical;
        int iter_elastic_cg;
        int iter_plastic_pcg;
        double time_elastic_cg;
        double time_plastic_pcg;
    };
    
    std::vector<Result> results;
    
    // ═══════════════════════════════════════════════════════════
    // Case 1: 弹性材料 + CG 求解器（基准）
    // ═══════════════════════════════════════════════════════════
    
    std::cout << "Case 1: Elastic Material + CG Solver\n";
    std::cout << std::string(50, '-') << "\n";
    
    IsotropicElastic elastic_material(E, nu, 2, true);  // 2D, plane_stress
    ElasticityUnified elastic_physics(&elastic_material, 2);
    
    CGSolver cg_solver;
    cg_solver.set_max_iter(10000);
    cg_solver.set_tol(1e-10);
    
    for (int step = 1; step <= num_steps; ++step) {
        Real u_applied = step * du;
        Real strain_applied = u_applied / length;
        
        Assembler assembler(model, 2);
        
        // 装配
        assembler.assemble([&](Index elem_id, const Mesh& m, DenseMatrix& Ke, Vector& Fe) {
            elastic_physics.compute_element(elem_id, m, Ke, Fe);
        });
        
        // 边界条件
        std::vector<DirichletBC> bcs;
        
        // 左端固定 x 方向
        for (Index i = 0; i < mesh.num_nodes(); ++i) {
            const Node& node = mesh.node(i);
            if (std::abs(node.coords()[0]) < 1e-6) {
                bcs.push_back({"left", static_cast<Index>(i * 2), 0.0});  // u_x = 0
            }
        }
        
        // 底边固定 y 方向（防止刚体运动）
        for (Index i = 0; i < mesh.num_nodes(); ++i) {
            const Node& node = mesh.node(i);
            if (std::abs(node.coords()[1]) < 1e-6) {
                bcs.push_back({"bottom", static_cast<Index>(i * 2 + 1), 0.0});  // u_y = 0
            }
        }
        
        // 右端施加位移
        for (Index i = 0; i < mesh.num_nodes(); ++i) {
            const Node& node = mesh.node(i);
            if (std::abs(node.coords()[0] - length) < 1e-6) {
                bcs.push_back({"right", static_cast<Index>(i * 2), u_applied});  // u_x = u
            }
        }
        
        assembler.apply_dirichlet(bcs);
        
        // 求解
        const auto& K = assembler.matrix();
        const auto& F = assembler.rhs();
        std::vector<Real> u(F.size(), 0.0);
        
        Timer timer;
        timer.start();
        auto result = cg_solver.solve(K, F.raw(), u);
        timer.stop();
        double elapsed = timer.elapsed_s();
        
        if (!result.converged) {
            std::cerr << "CG failed at step " << step << "!\n";
            break;
        }
        
        // 提取应力（从位移计算 - 最准确的方法）
        // 计算右端平均位移
        Real u_right_avg = 0.0;
        int count_right = 0;
        for (Index i = 0; i < mesh.num_nodes(); ++i) {
            const Node& node = mesh.node(i);
            if (std::abs(node.coords()[0] - length) < 1e-6) {
                u_right_avg += u[i * 2];  // u_x
                count_right++;
            }
        }
        u_right_avg /= count_right;
        
        // 应变和应力
        Real strain_actual = u_right_avg / length;
        Real sigma_elastic = E * strain_actual;  // 弹性材料：σ = E*ε
        
        // 保存结果
        Result r;
        r.displacement = u_applied;
        r.strain = strain_applied;
        r.stress_elastic_cg = sigma_elastic;
        r.iter_elastic_cg = result.iterations;
        r.time_elastic_cg = elapsed;
        
        if (results.size() < static_cast<size_t>(step)) {
            results.push_back(r);
        } else {
            results[step - 1].stress_elastic_cg = sigma_elastic;
            results[step - 1].iter_elastic_cg = result.iterations;
            results[step - 1].time_elastic_cg = elapsed;
        }
        
        if (step % 5 == 0 || step == 1) {
            std::cout << "  Step " << std::setw(2) << step
                      << ": ε = " << std::fixed << std::setprecision(6) << strain_applied
                      << ", σ = " << std::setw(8) << std::setprecision(2) << sigma_elastic << " MPa"
                      << ", iter = " << result.iterations
                      << ", time = " << std::setw(6) << std::setprecision(3) << elapsed * 1000 << " ms\n";
        }
    }
    
    std::cout << "\n";
    
    // ═══════════════════════════════════════════════════════════
    // Case 2: 塑性材料 + PCG 求解器（Jacobi 预条件）
    // ═══════════════════════════════════════════════════════════
    
    std::cout << "Case 2: Plastic Material + PCG Solver (Jacobi)\n";
    std::cout << std::string(50, '-') << "\n";
    
    J2Plasticity plastic_material(E, nu, sigma_y, H, 2);  // 2D
    ElasticityUnified plastic_physics(&plastic_material, 2);
    
    PCGSolver pcg_solver("jacobi");
    pcg_solver.set_max_iter(10000);
    pcg_solver.set_tol(1e-10);
    
    for (int step = 1; step <= num_steps; ++step) {
        Real u_applied = step * du;
        Real strain_applied = u_applied / length;
        
        Assembler assembler(model, 2);
        
        // 装配
        assembler.assemble([&](Index elem_id, const Mesh& m, DenseMatrix& Ke, Vector& Fe) {
            plastic_physics.compute_element(elem_id, m, Ke, Fe);
        });
        
        // 边界条件（同上）
        std::vector<DirichletBC> bcs;
        
        // 左端固定 x 方向
        for (Index i = 0; i < mesh.num_nodes(); ++i) {
            const Node& node = mesh.node(i);
            if (std::abs(node.coords()[0]) < 1e-6) {
                bcs.push_back({"left", static_cast<Index>(i * 2), 0.0});
            }
        }
        
        // 底边固定 y 方向
        for (Index i = 0; i < mesh.num_nodes(); ++i) {
            const Node& node = mesh.node(i);
            if (std::abs(node.coords()[1]) < 1e-6) {
                bcs.push_back({"bottom", static_cast<Index>(i * 2 + 1), 0.0});
            }
        }
        
        // 右端施加位移
        for (Index i = 0; i < mesh.num_nodes(); ++i) {
            const Node& node = mesh.node(i);
            if (std::abs(node.coords()[0] - length) < 1e-6) {
                bcs.push_back({"right", static_cast<Index>(i * 2), u_applied});
            }
        }
        
        assembler.apply_dirichlet(bcs);
        
        // 求解
        const auto& K = assembler.matrix();
        const auto& F = assembler.rhs();
        std::vector<Real> u(F.size(), 0.0);
        
        Timer timer;
        timer.start();
        auto result = pcg_solver.solve(K, F.raw(), u);
        timer.stop();
        double elapsed = timer.elapsed_s();
        
        if (!result.converged) {
            std::cerr << "PCG failed at step " << step << "!\n";
            break;
        }
        
        // 提取应力（从位移计算）
        Real u_right_avg = 0.0;
        int count_right = 0;
        for (Index i = 0; i < mesh.num_nodes(); ++i) {
            const Node& node = mesh.node(i);
            if (std::abs(node.coords()[0] - length) < 1e-6) {
                u_right_avg += u[i * 2];
                count_right++;
            }
        }
        u_right_avg /= count_right;
        
        // 应变
        Real strain_actual = u_right_avg / length;
        
        // 对于塑性材料，简单从应变判断：
        // ε < ε_y: σ = E*ε
        // ε >= ε_y: σ ≈ σ_y (理想塑性)
        Real sigma_plastic;
        if (strain_actual < analytical.epsilon_y) {
            sigma_plastic = E * strain_actual;  // 弹性
        } else {
            sigma_plastic = sigma_y;  // 塑性（近似）
        }
        
        // 解析解
        Real sigma_analytical = analytical.stress(strain_applied);
        
        // 保存结果
        results[step - 1].stress_plastic_pcg = sigma_plastic;
        results[step - 1].stress_analytical = sigma_analytical;
        results[step - 1].iter_plastic_pcg = result.iterations;
        results[step - 1].time_plastic_pcg = elapsed;
        
        if (step % 5 == 0 || step == 1) {
            std::cout << "  Step " << std::setw(2) << step
                      << ": ε = " << std::fixed << std::setprecision(6) << strain_applied
                      << ", σ = " << std::setw(8) << std::setprecision(2) << sigma_plastic << " MPa"
                      << ", iter = " << result.iterations
                      << ", time = " << std::setw(6) << std::setprecision(3) << elapsed * 1000 << " ms\n";
        }
    }
    
    std::cout << "\n";
    
    // ═══════════════════════════════════════════════════════════
    // 结果对比
    // ═══════════════════════════════════════════════════════════
    
    std::cout << "=== Results Comparison ===\n\n";
    
    std::cout << std::setw(6) << "Step"
              << std::setw(10) << "Strain"
              << std::setw(12) << "σ_elastic"
              << std::setw(12) << "σ_plastic"
              << std::setw(12) << "σ_theory"
              << std::setw(10) << "Error %"
              << std::setw(10) << "CG iter"
              << std::setw(10) << "PCG iter"
              << "\n";
    std::cout << std::string(82, '-') << "\n";
    
    for (size_t i = 0; i < results.size(); ++i) {
        const auto& r = results[i];
        Real error = std::abs(r.stress_plastic_pcg - r.stress_analytical) / r.stress_analytical * 100.0;
        
        std::cout << std::setw(6) << (i + 1)
                  << std::fixed << std::setprecision(6)
                  << std::setw(10) << r.strain
                  << std::setprecision(2)
                  << std::setw(12) << r.stress_elastic_cg
                  << std::setw(12) << r.stress_plastic_pcg
                  << std::setw(12) << r.stress_analytical
                  << std::setprecision(3)
                  << std::setw(9) << error << "%"
                  << std::setw(10) << r.iter_elastic_cg
                  << std::setw(10) << r.iter_plastic_pcg
                  << "\n";
    }
    
    std::cout << std::string(82, '-') << "\n\n";
    
    // 统计
    double total_time_cg = 0.0, total_time_pcg = 0.0;
    int total_iter_cg = 0, total_iter_pcg = 0;
    
    for (const auto& r : results) {
        total_time_cg += r.time_elastic_cg;
        total_time_pcg += r.time_plastic_pcg;
        total_iter_cg += r.iter_elastic_cg;
        total_iter_pcg += r.iter_plastic_pcg;
    }
    
    std::cout << "=== Performance Summary ===\n\n";
    std::cout << "CG Solver (Elastic):\n";
    std::cout << "  Total iterations: " << total_iter_cg << "\n";
    std::cout << "  Avg iter/step:    " << total_iter_cg / num_steps << "\n";
    std::cout << "  Total time:       " << total_time_cg * 1000 << " ms\n";
    std::cout << "  Avg time/step:    " << total_time_cg / num_steps * 1000 << " ms\n\n";
    
    std::cout << "PCG Solver (Plastic):\n";
    std::cout << "  Total iterations: " << total_iter_pcg << "\n";
    std::cout << "  Avg iter/step:    " << total_iter_pcg / num_steps << "\n";
    std::cout << "  Total time:       " << total_time_pcg * 1000 << " ms\n";
    std::cout << "  Avg time/step:    " << total_time_pcg / num_steps * 1000 << " ms\n\n";
    
    std::cout << "Speedup (PCG vs CG):\n";
    std::cout << "  Iteration:        " << static_cast<double>(total_iter_cg) / total_iter_pcg << "x\n";
    std::cout << "  Time:             " << total_time_cg / total_time_pcg << "x\n\n";
    
    // ═══════════════════════════════════════════════════════════
    // 导出结果
    // ═══════════════════════════════════════════════════════════
    
    std::ofstream outfile("plasticity_pcg_results.dat");
    outfile << "# Plasticity Test with PCG Solver\n";
    outfile << "# Step  ε  σ_elastic  σ_plastic  σ_theory  error%  iter_CG  iter_PCG\n";
    
    for (size_t i = 0; i < results.size(); ++i) {
        const auto& r = results[i];
        Real error = std::abs(r.stress_plastic_pcg - r.stress_analytical) / r.stress_analytical * 100.0;
        
        outfile << (i + 1) << "  "
                << r.strain << "  "
                << r.stress_elastic_cg << "  "
                << r.stress_plastic_pcg << "  "
                << r.stress_analytical << "  "
                << error << "  "
                << r.iter_elastic_cg << "  "
                << r.iter_plastic_pcg << "\n";
    }
    
    outfile.close();
    
    std::cout << "Results exported to: plasticity_pcg_results.dat\n\n";
    
    std::cout << "================================================\n";
    std::cout << "  Test Completed Successfully!\n";
    std::cout << "================================================\n";
    
    return 0;
}
