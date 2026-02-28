/**
 * test_plasticity_verification.cpp - 塑性材料验证测试
 * 
 * 简化测试：验证 J2 塑性材料的基本行为
 * - 单元测试级别
 * - 与解析解对比
 */

#include "material/j2_plasticity.h"
#include "material/isotropic_elastic.h"
#include "math/dense_matrix.h"
#include "core/logger.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace fem;
using namespace fem::constitutive;

int main() {
    Logger::instance().set_level(LogLevel::WARN);
    
    std::cout << "======================================\n";
    std::cout << "  J2 Plasticity Verification Test\n";
    std::cout << "======================================\n\n";
    
    // ═══════════════════════════════════════════════════════════
    // 材料参数
    // ═══════════════════════════════════════════════════════════
    
    Real E = 200e3;        // MPa
    Real nu = 0.3;
    Real sigma_y = 250.0;  // MPa
    Real H = 0.0;          // 理想塑性
    
    std::cout << "=== Material Parameters ===\n";
    std::cout << "E   = " << E << " MPa\n";
    std::cout << "ν   = " << nu << "\n";
    std::cout << "σ_y = " << sigma_y << " MPa\n";
    std::cout << "H   = " << H << " MPa (perfect plastic)\n\n";
    
    // 创建材料
    J2Plasticity plastic(E, nu, sigma_y, H, 2);  // 2D, plane stress
    IsotropicElastic elastic(E, nu, 2, true);     // 2D, plane stress
    
    // ═══════════════════════════════════════════════════════════
    // Test 1: 单轴拉伸（解析解）
    // ═══════════════════════════════════════════════════════════
    
    std::cout << "=== Test 1: Uniaxial Tension ===\n\n";
    
    Real epsilon_y = sigma_y / E;  // 屈服应变
    
    std::cout << "Yield strain: ε_y = " << epsilon_y << "\n";
    std::cout << "Analytical solution:\n";
    std::cout << "  ε < ε_y:  σ = E*ε\n";
    std::cout << "  ε > ε_y:  σ = σ_y (H=0)\n\n";
    
    std::cout << std::setw(10) << "ε"
              << std::setw(12) << "σ_elastic"
              << std::setw(12) << "σ_plastic"
              << std::setw(12) << "σ_theory"
              << std::setw(12) << "Error %\n";
    std::cout << std::string(58, '-') << "\n";
    
    std::vector<Real> strain_list, stress_elastic_list, stress_plastic_list;
    
    int n_steps = 20;
    Real max_strain = 0.003;  // 3 * ε_y
    
    // 塑性材料需要增量加载
    Vector stress_plastic_current(3, 0.0);
    StateVariables state_plastic = plastic.createState();  // 初始化状态变量
    
    for (int i = 1; i <= n_steps; ++i) {
        Real epsilon = (max_strain / n_steps) * i;
        
        // 应变增量
        Real d_epsilon = max_strain / n_steps;
        Vector strain_inc(3);
        strain_inc[0] = d_epsilon;  // Δε_xx
        strain_inc[1] = 0.0;        // Δε_yy
        strain_inc[2] = 0.0;        // Δγ_xy
        
        // 总应变 (用于弹性材料)
        Vector strain_total(3);
        strain_total[0] = epsilon;
        strain_total[1] = 0.0;
        strain_total[2] = 0.0;
        
        // --- Elastic material (直接计算) ---
        DenseMatrix D_elastic;
        StateVariables state_elastic;
        elastic.computeTangent(strain_total, D_elastic, state_elastic);
        
        Vector stress_elastic = D_elastic * strain_total;
        Real sigma_elastic = stress_elastic[0];  // σ_xx
        
        // --- Plastic material (增量更新) ---
        plastic.computeStress(strain_inc, stress_plastic_current, state_plastic);
        
        Real sigma_plastic = stress_plastic_current[0];  // σ_xx
        
        // --- Analytical solution ---
        Real sigma_theory = (epsilon < epsilon_y) ? E * epsilon : sigma_y;
        
        // 相对误差
        Real error = std::abs(sigma_plastic - sigma_theory) / sigma_theory * 100.0;
        
        std::cout << std::fixed << std::setprecision(6)
                  << std::setw(10) << epsilon
                  << std::setw(12) << sigma_elastic
                  << std::setw(12) << sigma_plastic
                  << std::setw(12) << sigma_theory
                  << std::setw(11) << error << "%\n";
        
        strain_list.push_back(epsilon);
        stress_elastic_list.push_back(sigma_elastic);
        stress_plastic_list.push_back(sigma_plastic);
    }
    
    std::cout << std::string(58, '-') << "\n\n";
    
    // ═══════════════════════════════════════════════════════════
    // Test 2: 屈服判据验证
    // ═══════════════════════════════════════════════════════════
    
    std::cout << "=== Test 2: Yield Criterion ===\n\n";
    
    std::cout << "von Mises yield function: f = √(3*J2) - σ_y\n";
    std::cout << "  J2 = 1/2 * s:s (s = deviatoric stress)\n";
    std::cout << "  For uniaxial tension: f = |σ| - σ_y\n\n";
    
    // 测试屈服面
    for (Real factor : {0.5, 0.9, 1.0, 1.1, 1.5}) {
        Real epsilon = factor * epsilon_y;
        
        Vector strain(3);
        strain[0] = epsilon;
        strain[1] = 0.0;
        strain[2] = 0.0;
        
        DenseMatrix D;
        StateVariables state;
        plastic.computeTangent(strain, D, state);
        
        Vector stress = D * strain;
        Real sigma = stress[0];
        
        // von Mises 等效应力
        Real sigma_vm = std::abs(sigma);  // 单轴拉伸
        
        // 屈服函数值
        Real f = sigma_vm - sigma_y;
        
        std::cout << "ε/ε_y = " << std::setw(4) << factor
                  << "  →  σ = " << std::setw(8) << sigma << " MPa"
                  << ",  σ_vm = " << std::setw(8) << sigma_vm << " MPa"
                  << ",  f = " << std::setw(8) << f << " MPa";
        
        if (f < -1e-3) {
            std::cout << "  [Elastic]\n";
        } else if (std::abs(f) < 1e-3) {
            std::cout << "  [On yield surface]\n";
        } else {
            std::cout << "  [Plastic]\n";
        }
    }
    
    std::cout << "\n";
    
    // ═══════════════════════════════════════════════════════════
    // 导出结果
    // ═══════════════════════════════════════════════════════════
    
    std::ofstream outfile("plasticity_verification.dat");
    outfile << "# J2 Plasticity Verification\n";
    outfile << "#   E = " << E << " MPa\n";
    outfile << "#  nu = " << nu << "\n";
    outfile << "# σ_y = " << sigma_y << " MPa\n";
    outfile << "#   H = " << H << " MPa\n";
    outfile << "# ε_y = " << epsilon_y << "\n";
    outfile << "#\n";
    outfile << "# ε           σ_elastic    σ_plastic    σ_theory\n";
    
    for (size_t i = 0; i < strain_list.size(); ++i) {
        Real sigma_theory = (strain_list[i] < epsilon_y) ? 
                            E * strain_list[i] : sigma_y;
        
        outfile << strain_list[i] << "  "
                << stress_elastic_list[i] << "  "
                << stress_plastic_list[i] << "  "
                << sigma_theory << "\n";
    }
    
    outfile.close();
    
    std::cout << "Results exported to: plasticity_verification.dat\n\n";
    
    // ═══════════════════════════════════════════════════════════
    // 总结
    // ═══════════════════════════════════════════════════════════
    
    std::cout << "=== Summary ===\n\n";
    
    // 找到屈服点
    size_t yield_idx = 0;
    for (size_t i = 0; i < strain_list.size(); ++i) {
        if (std::abs(stress_plastic_list[i] - sigma_y) / sigma_y < 0.01) {
            yield_idx = i;
            break;
        }
    }
    
    if (yield_idx > 0) {
        std::cout << "Yielding observed at:\n";
        std::cout << "  ε_FEM   = " << strain_list[yield_idx] << "\n";
        std::cout << "  ε_exact = " << epsilon_y << "\n";
        std::cout << "  Error   = " << std::abs(strain_list[yield_idx] - epsilon_y) / epsilon_y * 100 << "%\n\n";
    }
    
    // 最终应力
    std::cout << "At maximum strain (ε = " << max_strain << "):\n";
    std::cout << "  σ_elastic = " << stress_elastic_list.back() << " MPa\n";
    std::cout << "  σ_plastic = " << stress_plastic_list.back() << " MPa\n";
    std::cout << "  σ_theory  = " << sigma_y << " MPa (H=0)\n";
    std::cout << "  Error     = " << std::abs(stress_plastic_list.back() - sigma_y) / sigma_y * 100 << "%\n\n";
    
    std::cout << "======================================\n";
    std::cout << "  Verification Completed!\n";
    std::cout << "======================================\n";
    
    return 0;
}
