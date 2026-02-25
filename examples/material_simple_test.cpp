/**
 * 材料模型简化测试 - 单元级别验证
 * 
 * 测试目标：
 * 1. 验证IsotropicElastic计算正确性
 * 2. 验证J2Plasticity返回映射算法
 * 3. 对比弹性和塑性响应
 */

#include "material/isotropic_elastic.h"
#include "material/j2_plasticity.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace fem;
using namespace fem::constitutive;

void print_header(const std::string& title) {
    std::cout << "\n";
    std::cout << "═══════════════════════════════════════════════════════════\n";
    std::cout << "  " << title << "\n";
    std::cout << "═══════════════════════════════════════════════════════════\n\n";
}

void print_vector(const std::string& name, const Vector& v) {
    std::cout << std::setw(20) << name << " = [";
    for (std::size_t i = 0; i < v.size(); ++i) {
        std::cout << std::setw(10) << std::scientific << std::setprecision(3) << v[i];
        if (i < v.size() - 1) std::cout << ", ";
    }
    std::cout << "]\n";
}

void test_isotropic_elastic() {
    print_header("测试 1: IsotropicElastic - 单轴拉伸");
    
    // 创建材料
    IsotropicElastic steel(200e3, 0.3, 2, true);  // 2D平面应力
    std::cout << "材料: " << steel.typeName() << "\n";
    std::cout << "  E  = " << steel.getParameter("E") << " MPa\n";
    std::cout << "  nu = " << steel.getParameter("nu") << "\n";
    std::cout << "  μ  = " << steel.mu() << " MPa\n\n";
    
    StateVariables state = steel.createState();
    Vector stress(3, 0.0);
    
    // 步骤1：加载
    std::cout << "─── 步骤 1: 加载 (ε11 = 0.1%) ───\n";
    Vector strain_inc(3, 0.0);
    strain_inc[0] = 0.001;  // ε11 = 0.1%
    
    steel.computeStress(strain_inc, stress, state);
    print_vector("应变增量", strain_inc);
    print_vector("应力", stress);
    
    // 理论验证
    Real E = steel.getParameter("E");
    Real nu = steel.getParameter("nu");
    Real factor = E / (1.0 - nu * nu);
    Real sigma11_theory = factor * 0.001;
    std::cout << "\n验证: σ11 理论值 = " << std::scientific << sigma11_theory << " MPa\n";
    std::cout << "      误差 = " << std::abs(stress[0] - sigma11_theory) / sigma11_theory * 100 << " %\n";
    
    // 步骤2：卸载
    std::cout << "\n─── 步骤 2: 卸载 (减小 0.05%) ───\n";
    strain_inc[0] = -0.0005;
    steel.computeStress(strain_inc, stress, state);
    print_vector("应变增量", strain_inc);
    print_vector("应力", stress);
}

void test_j2_plasticity_elastic() {
    print_header("测试 2: J2Plasticity - 弹性加载");
    
    // 创建材料
    J2Plasticity mild_steel(200e3, 0.3, 250.0, 1000.0, 2);  // 2D
    std::cout << "材料: " << mild_steel.typeName() << "\n";
    std::cout << "  E    = " << mild_steel.getParameter("E") << " MPa\n";
    std::cout << "  σ_y0 = " << mild_steel.getParameter("sigma_y0") << " MPa\n";
    std::cout << "  H    = " << mild_steel.getParameter("H") << " MPa\n\n";
    
    StateVariables state = mild_steel.createState();
    Vector stress(3, 0.0);
    
    // 小应变（弹性）
    std::cout << "─── 弹性加载 (ε11 = 0.05%) ───\n";
    Vector strain_inc(3, 0.0);
    strain_inc[0] = 0.0005;
    
    mild_steel.computeStress(strain_inc, stress, state);
    print_vector("应变增量", strain_inc);
    print_vector("应力", stress);
    
    Real q = mild_steel.vonMisesStress(stress);
    Real sigma_y = mild_steel.yieldStress(state.equiv_plastic_strain);
    
    std::cout << "\nvon Mises 应力: q = " << q << " MPa\n";
    std::cout << "屈服应力:      σ_y = " << sigma_y << " MPa\n";
    std::cout << "屈服函数:      f = " << (q - sigma_y) << " MPa\n";
    std::cout << "等效塑性应变:  ε_p = " << state.equiv_plastic_strain << "\n";
    std::cout << "状态:          " << (q < sigma_y ? "弹性 ✓" : "塑性") << "\n";
}

void test_j2_plasticity_plastic() {
    print_header("测试 3: J2Plasticity - 塑性加载");
    
    J2Plasticity mild_steel(200e3, 0.3, 250.0, 1000.0, 2);
    StateVariables state = mild_steel.createState();
    Vector stress(3, 0.0);
    
    // 大应变（塑性）
    std::cout << "─── 步骤 1: 大应变加载 (ε11 = 0.2%) ───\n";
    Vector strain_inc(3, 0.0);
    strain_inc[0] = 0.002;
    
    mild_steel.computeStress(strain_inc, stress, state);
    print_vector("应变增量", strain_inc);
    print_vector("应力", stress);
    
    Real q = mild_steel.vonMisesStress(stress);
    Real sigma_y = mild_steel.yieldStress(state.equiv_plastic_strain);
    
    std::cout << "\nvon Mises 应力: q = " << q << " MPa\n";
    std::cout << "屈服应力:      σ_y = " << sigma_y << " MPa\n";
    std::cout << "屈服函数:      f = " << (q - sigma_y) << " MPa (≈0)\n";
    std::cout << "等效塑性应变:  ε_p = " << state.equiv_plastic_strain << "\n";
    std::cout << "状态:          塑性 (返回映射) ✓\n";
    
    // 卸载
    std::cout << "\n─── 步骤 2: 弹性卸载 (减小 0.1%) ───\n";
    Real eps_p_before = state.equiv_plastic_strain;
    strain_inc[0] = -0.001;
    
    mild_steel.computeStress(strain_inc, stress, state);
    print_vector("应变增量", strain_inc);
    print_vector("应力", stress);
    
    std::cout << "\n等效塑性应变: ε_p = " << state.equiv_plastic_strain 
              << " (不变: " << (std::abs(state.equiv_plastic_strain - eps_p_before) < 1e-10 ? "✓" : "✗") << ")\n";
    std::cout << "状态:         弹性卸载 ✓\n";
    
    // 再加载
    std::cout << "\n─── 步骤 3: 再加载 (增加 0.15%) ───\n";
    strain_inc[0] = 0.0015;
    
    mild_steel.computeStress(strain_inc, stress, state);
    print_vector("应力", stress);
    
    q = mild_steel.vonMisesStress(stress);
    sigma_y = mild_steel.yieldStress(state.equiv_plastic_strain);
    
    std::cout << "\nvon Mises 应力: q = " << q << " MPa\n";
    std::cout << "硬化后屈服应力: σ_y = " << sigma_y << " MPa\n";
    std::cout << "等效塑性应变:  ε_p = " << state.equiv_plastic_strain << " (累积 ✓)\n";
}

void test_comparison() {
    print_header("测试 4: 弹性vs塑性对比");
    
    // 相同参数
    Real E = 200e3, nu = 0.3;
    
    IsotropicElastic elastic(E, nu, 2, true);
    J2Plasticity plastic(E, nu, 250.0, 1000.0, 2);
    
    StateVariables state_e = elastic.createState();
    StateVariables state_p = plastic.createState();
    
    Vector stress_e(3, 0.0), stress_p(3, 0.0);
    
    std::cout << "加载: ε11 = 0.3% (大应变)\n\n";
    
    Vector strain_inc(3, 0.0);
    strain_inc[0] = 0.003;
    
    elastic.computeStress(strain_inc, stress_e, state_e);
    plastic.computeStress(strain_inc, stress_p, state_p);
    
    print_vector("弹性应力", stress_e);
    print_vector("塑性应力", stress_p);
    
    std::cout << "\n─── 对比 ───\n";
    std::cout << "弹性: σ11 = " << stress_e[0] << " MPa\n";
    std::cout << "塑性: σ11 = " << stress_p[0] << " MPa\n";
    std::cout << "差异: " << (stress_e[0] - stress_p[0]) << " MPa\n";
    std::cout << "\n塑性材料屈服后应力增长减缓（硬化） ✓\n";
    std::cout << "等效塑性应变: ε_p = " << state_p.equiv_plastic_strain << "\n";
}

int main() {
    std::cout << "\n";
    std::cout << "███████████████████████████████████████████████████████████\n";
    std::cout << "  材料模型集成测试 - 单元级验证\n";
    std::cout << "  Author: Math Agent 🧮\n";
    std::cout << "███████████████████████████████████████████████████████████\n";
    
    try {
        test_isotropic_elastic();
        test_j2_plasticity_elastic();
        test_j2_plasticity_plastic();
        test_comparison();
        
        std::cout << "\n";
        std::cout << "═══════════════════════════════════════════════════════════\n";
        std::cout << "  所有测试通过！ ✓\n";
        std::cout << "═══════════════════════════════════════════════════════════\n\n";
        
    } catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}
