/**
 * 材料本构模型演示程序
 * 
 * 演示 IsotropicElastic 和 J2Plasticity 的使用
 */

#include "material/isotropic_elastic.h"
#include "material/j2_plasticity.h"
#include <iostream>
#include <iomanip>

using namespace fem;
using namespace fem::constitutive;

void printVector(const std::string& name, const Vector& vec) {
    std::cout << std::setw(20) << name << " = [";
    for (std::size_t i = 0; i < vec.size(); ++i) {
        std::cout << std::setw(10) << std::scientific << std::setprecision(3) << vec[i];
        if (i < vec.size() - 1) std::cout << ", ";
    }
    std::cout << "]\n";
}

void demoIsotropicElastic() {
    std::cout << "\n";
    std::cout << "═══════════════════════════════════════════════════════════\n";
    std::cout << "  IsotropicElastic - 各向同性弹性材料演示\n";
    std::cout << "═══════════════════════════════════════════════════════════\n\n";
    
    // 创建材料：钢材 (E=200GPa, nu=0.3)
    IsotropicElastic steel(200e3, 0.3, 3);
    std::cout << "材料: " << steel.typeName() << "\n";
    std::cout << "  E  = " << steel.getParameter("E") << " MPa\n";
    std::cout << "  nu = " << steel.getParameter("nu") << "\n";
    std::cout << "  λ  = " << steel.lambda() << " MPa\n";
    std::cout << "  μ  = " << steel.mu() << " MPa\n\n";
    
    // 初始化状态
    StateVariables state = steel.createState();
    Vector stress(6, 0.0);
    
    // 单轴拉伸试验
    std::cout << "─── 单轴拉伸 (ε11 = 0.1%) ───\n";
    Vector strain_inc(6, 0.0);
    strain_inc[0] = 0.001;  // ε11 = 0.1%
    
    steel.computeStress(strain_inc, stress, state);
    printVector("应变增量", strain_inc);
    printVector("应力", stress);
    
    Vector total_strain = strain_inc;
    Real energy = steel.strainEnergy(total_strain, state);
    std::cout << std::setw(20) << "应变能密度" << " = " 
              << std::scientific << energy << " MPa\n";
    
    // 继续加载
    std::cout << "\n─── 继续加载 (再增加 0.1%) ───\n";
    steel.computeStress(strain_inc, stress, state);
    printVector("应力", stress);
    
    // 纯剪切
    std::cout << "\n─── 纯剪切 (γ12 = 0.2%) ───\n";
    stress.zero();
    strain_inc.zero();
    strain_inc[3] = 0.002;  // γ12 = 0.2%
    
    steel.computeStress(strain_inc, stress, state);
    printVector("应变增量", strain_inc);
    printVector("应力", stress);
    std::cout << "  验证: σ12 = G * γ12 = " << steel.mu() << " * 0.002 = " 
              << steel.mu() * 0.002 << " MPa\n";
}

void demoJ2Plasticity() {
    std::cout << "\n";
    std::cout << "═══════════════════════════════════════════════════════════\n";
    std::cout << "  J2Plasticity - von Mises 塑性材料演示\n";
    std::cout << "═══════════════════════════════════════════════════════════\n\n";
    
    // 创建材料：软钢 (σy0=250MPa, H=1000MPa)
    J2Plasticity mild_steel(200e3, 0.3, 250.0, 1000.0, 3);
    std::cout << "材料: " << mild_steel.typeName() << "\n";
    std::cout << "  E      = " << mild_steel.getParameter("E") << " MPa\n";
    std::cout << "  nu     = " << mild_steel.getParameter("nu") << "\n";
    std::cout << "  σ_y0   = " << mild_steel.getParameter("sigma_y0") << " MPa\n";
    std::cout << "  H      = " << mild_steel.getParameter("H") << " MPa\n\n";
    
    StateVariables state = mild_steel.createState();
    Vector stress(6, 0.0);
    
    // 第1步：弹性加载
    std::cout << "─── Step 1: 弹性加载 (ε11 = 0.05%) ───\n";
    Vector strain_inc(6, 0.0);
    strain_inc[0] = 0.0005;
    
    mild_steel.computeStress(strain_inc, stress, state);
    Real q = mild_steel.vonMisesStress(stress);
    Real sigma_y = mild_steel.yieldStress(state.equiv_plastic_strain);
    
    printVector("应变增量", strain_inc);
    printVector("应力", stress);
    std::cout << std::setw(20) << "von Mises 应力" << " = " << q << " MPa\n";
    std::cout << std::setw(20) << "当前屈服应力" << " = " << sigma_y << " MPa\n";
    std::cout << std::setw(20) << "等效塑性应变" << " = " << state.equiv_plastic_strain << "\n";
    std::cout << std::setw(20) << "状态" << " = " << (q < sigma_y ? "弹性" : "塑性") << "\n";
    
    // 第2步：塑性加载
    std::cout << "\n─── Step 2: 塑性加载 (再增加 0.15%) ───\n";
    strain_inc[0] = 0.0015;
    
    mild_steel.computeStress(strain_inc, stress, state);
    q = mild_steel.vonMisesStress(stress);
    sigma_y = mild_steel.yieldStress(state.equiv_plastic_strain);
    
    printVector("应变增量", strain_inc);
    printVector("应力", stress);
    std::cout << std::setw(20) << "von Mises 应力" << " = " << q << " MPa\n";
    std::cout << std::setw(20) << "当前屈服应力" << " = " << sigma_y << " MPa\n";
    std::cout << std::setw(20) << "等效塑性应变" << " = " << state.equiv_plastic_strain << "\n";
    std::cout << std::setw(20) << "屈服函数" << " = " << (q - sigma_y) << " MPa (≈0)\n";
    std::cout << std::setw(20) << "状态" << " = 塑性（返回映射）\n";
    
    // 第3步：弹性卸载
    std::cout << "\n─── Step 3: 弹性卸载 (减小 0.1%) ───\n";
    strain_inc[0] = -0.001;
    
    Real eps_p_before = state.equiv_plastic_strain;
    mild_steel.computeStress(strain_inc, stress, state);
    q = mild_steel.vonMisesStress(stress);
    
    printVector("应变增量", strain_inc);
    printVector("应力", stress);
    std::cout << std::setw(20) << "von Mises 应力" << " = " << q << " MPa\n";
    std::cout << std::setw(20) << "等效塑性应变" << " = " << state.equiv_plastic_strain << " (不变)\n";
    std::cout << std::setw(20) << "状态" << " = 弹性卸载\n";
    
    // 第4步：再加载
    std::cout << "\n─── Step 4: 再加载 (增加 0.15%) ───\n";
    strain_inc[0] = 0.0015;
    
    mild_steel.computeStress(strain_inc, stress, state);
    q = mild_steel.vonMisesStress(stress);
    sigma_y = mild_steel.yieldStress(state.equiv_plastic_strain);
    
    printVector("应力", stress);
    std::cout << std::setw(20) << "von Mises 应力" << " = " << q << " MPa\n";
    std::cout << std::setw(20) << "等效塑性应变" << " = " << state.equiv_plastic_strain 
              << " (累积)\n";
    std::cout << std::setw(20) << "硬化后屈服应力" << " = " << sigma_y << " MPa\n";
}

int main() {
    std::cout << "\n";
    std::cout << "███████████████████████████████████████████████████████████\n";
    std::cout << "  材料本构模型演示程序\n";
    std::cout << "  Author: Math Agent 🧮\n";
    std::cout << "███████████████████████████████████████████████████████████\n";
    
    demoIsotropicElastic();
    demoJ2Plasticity();
    
    std::cout << "\n";
    std::cout << "═══════════════════════════════════════════════════════════\n";
    std::cout << "  演示完成！\n";
    std::cout << "═══════════════════════════════════════════════════════════\n\n";
    
    return 0;
}
