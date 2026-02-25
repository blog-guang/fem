#include "j2_plasticity.h"
#include <cmath>
#include <stdexcept>
#include <iostream>

namespace fem {
namespace constitutive {

J2Plasticity::J2Plasticity(Real E, Real nu, 
                           Real sigma_y0, 
                           Real H,
                           int dimension)
    : Material(dimension == 3 ? 6 : 3),
      dimension_(dimension)
{
    setParameter("E", E);
    setParameter("nu", nu);
    setParameter("sigma_y0", sigma_y0);
    setParameter("H", H);
    validateParameters();
}

void J2Plasticity::computeStress(
    const Vector& strain_inc, 
    Vector& stress, 
    StateVariables& state
) {
    // 初始化应力
    if (stress.size() != strain_size_) {
        stress.resize(strain_size_, 0.0);
    }
    
    // 1. 弹性预测：σ_trial = σ_n + D^e : Δε
    DenseMatrix D_elastic = elasticTensor();
    Vector delta_stress(strain_size_, 0.0);
    
    for (std::size_t i = 0; i < strain_size_; ++i) {
        for (std::size_t j = 0; j < strain_size_; ++j) {
            delta_stress[i] += D_elastic(i, j) * strain_inc[j];
        }
    }
    
    Vector stress_trial = stress + delta_stress;
    
    // 2. 检查屈服
    Real q_trial = vonMisesStress(stress_trial);
    Real sigma_y = yieldStress(state.equiv_plastic_strain);
    Real f_trial = q_trial - sigma_y;
    
    if (f_trial <= 0.0) {
        // 弹性步：应力在屈服面内
        stress = stress_trial;
    } else {
        // 塑性步：返回映射
        DenseMatrix D_ep;  // 弹塑性切线（本函数不需要，但returnMapping需要）
        returnMapping(stress_trial, state, D_ep);
        stress = stress_trial;  // stress_trial已被修改
    }
}

void J2Plasticity::computeTangent(
    const Vector& /*strain*/,
    DenseMatrix& D_mat,
    const StateVariables& state
) {
    // 简化实现：返回弹性刚度
    // 完整实现需要判断当前是否处于塑性状态并返回一致性切线
    // TODO: 基于当前应力状态判断是否需要弹塑性切线
    
    D_mat = elasticTensor();
    
    // 注意：更准确的实现需要存储上一步的塑性乘子
    // 这里作为简化，返回弹性刚度（适用于显式算法）
}

Real J2Plasticity::strainEnergy(
    const Vector& strain,
    const StateVariables& state
) const {
    // 总应变能 = 弹性应变能 + 塑性耗散
    // 这里仅计算弹性部分：Ψ_e = 0.5 * ε_e : D : ε_e
    // ε_e = ε - ε_p
    
    Vector elastic_strain = strain;
    for (std::size_t i = 0; i < strain_size_; ++i) {
        elastic_strain[i] -= state.plastic_strain[i];
    }
    
    DenseMatrix D = elasticTensor();
    Real energy = 0.0;
    for (std::size_t i = 0; i < strain_size_; ++i) {
        for (std::size_t j = 0; j < strain_size_; ++j) {
            energy += 0.5 * elastic_strain[i] * D(i, j) * elastic_strain[j];
        }
    }
    return energy;
}

StateVariables J2Plasticity::createState() const {
    StateVariables state(strain_size_);
    // 塑性应变和等效塑性应变初始化为0（已在构造函数中完成）
    return state;
}

std::string J2Plasticity::typeName() const {
    return dimension_ == 3 ? "J2Plasticity (3D)" : "J2Plasticity (2D)";
}

void J2Plasticity::validateParameters() const {
    Real E = getParameter("E");
    Real nu = getParameter("nu");
    Real sigma_y0 = getParameter("sigma_y0");
    Real H = getParameter("H");
    
    if (E <= 0.0) {
        throw std::invalid_argument("Young's modulus E must be positive");
    }
    if (nu <= -1.0 || nu >= 0.5) {
        throw std::invalid_argument("Poisson's ratio must be in (-1, 0.5)");
    }
    if (sigma_y0 <= 0.0) {
        throw std::invalid_argument("Initial yield stress sigma_y0 must be positive");
    }
    if (H < 0.0) {
        throw std::invalid_argument("Hardening modulus H must be non-negative");
    }
}

Real J2Plasticity::vonMisesStress(const Vector& stress) const {
    Vector s = deviatoricStress(stress);
    
    // q = √(3/2 * s:s)
    Real s_norm_sq = 0.0;
    if (dimension_ == 3) {
        // 3D: s:s = s11^2 + s22^2 + s33^2 + 2*(s12^2 + s23^2 + s13^2)
        s_norm_sq = s[0]*s[0] + s[1]*s[1] + s[2]*s[2] 
                  + 2.0 * (s[3]*s[3] + s[4]*s[4] + s[5]*s[5]);
    } else {
        // 2D: s:s = s11^2 + s22^2 + 2*s12^2
        s_norm_sq = s[0]*s[0] + s[1]*s[1] + 2.0 * s[2]*s[2];
    }
    
    return std::sqrt(1.5 * s_norm_sq);
}

Real J2Plasticity::yieldFunction(Real equiv_stress, Real equiv_plastic_strain) const {
    return equiv_stress - yieldStress(equiv_plastic_strain);
}

Real J2Plasticity::yieldStress(Real equiv_plastic_strain) const {
    Real sigma_y0 = getParameter("sigma_y0");
    Real H = getParameter("H");
    return sigma_y0 + H * equiv_plastic_strain;
}

void J2Plasticity::returnMapping(
    Vector& stress_trial,
    StateVariables& state,
    DenseMatrix& D_ep
) {
    // 径向返回映射算法（Radial Return Mapping）
    
    Real E = getParameter("E");
    Real nu = getParameter("nu");
    Real G = E / (2.0 * (1.0 + nu));  // 剪切模量
    Real H = getParameter("H");
    
    // 1. 分解试探应力
    Real p = hydrostaticPressure(stress_trial);
    Vector s_trial = deviatoricStress(stress_trial);
    Real q_trial = vonMisesStress(stress_trial);
    
    // 2. 计算塑性乘子增量（Newton迭代）
    Real sigma_y = yieldStress(state.equiv_plastic_strain);
    Real delta_gamma = (q_trial - sigma_y) / (3.0 * G + H);
    
    if (delta_gamma < 0.0) {
        delta_gamma = 0.0;  // 安全检查
    }
    
    // 3. 修正应力：径向返回
    Real factor = 1.0 - (3.0 * G * delta_gamma) / q_trial;
    if (factor < 0.0) factor = 0.0;  // 防止数值问题
    
    Vector s_corrected(strain_size_);
    for (std::size_t i = 0; i < strain_size_; ++i) {
        s_corrected[i] = factor * s_trial[i];
    }
    
    // 4. 重构应力：σ = s + p*I
    if (dimension_ == 3) {
        stress_trial[0] = s_corrected[0] + p;
        stress_trial[1] = s_corrected[1] + p;
        stress_trial[2] = s_corrected[2] + p;
        stress_trial[3] = s_corrected[3];
        stress_trial[4] = s_corrected[4];
        stress_trial[5] = s_corrected[5];
    } else {
        stress_trial[0] = s_corrected[0] + p;
        stress_trial[1] = s_corrected[1] + p;
        stress_trial[2] = s_corrected[2];
    }
    
    // 5. 更新状态变量
    // 塑性应变增量：Δε_p = √(3/2) * Δγ * n
    // 其中 n = s / ||s|| 是流动方向
    Real sqrt_3_2 = std::sqrt(1.5);
    Real q_new = vonMisesStress(stress_trial);
    Real norm_s = q_new / sqrt_3_2;
    
    if (norm_s > 1e-10) {
        for (std::size_t i = 0; i < strain_size_; ++i) {
            Real flow_dir = s_corrected[i] / norm_s;
            state.plastic_strain[i] += sqrt_3_2 * delta_gamma * flow_dir;
        }
    }
    
    // 更新等效塑性应变
    state.equiv_plastic_strain += delta_gamma;
    
    // 6. 计算一致性切线刚度
    D_ep = consistentTangent(stress_trial, delta_gamma, state.equiv_plastic_strain);
}

Vector J2Plasticity::deviatoricStress(const Vector& stress) const {
    Real p = hydrostaticPressure(stress);
    
    Vector s = stress;
    if (dimension_ == 3) {
        s[0] -= p;
        s[1] -= p;
        s[2] -= p;
        // 剪应力分量保持不变
    } else {
        s[0] -= p;
        s[1] -= p;
        // s[2] (剪应力) 保持不变
    }
    
    return s;
}

Real J2Plasticity::hydrostaticPressure(const Vector& stress) const {
    if (dimension_ == 3) {
        return (stress[0] + stress[1] + stress[2]) / 3.0;
    } else {
        return (stress[0] + stress[1]) / 2.0;  // 2D近似
    }
}

DenseMatrix J2Plasticity::elasticTensor() const {
    Real E = getParameter("E");
    Real nu = getParameter("nu");
    
    if (dimension_ == 3) {
        // 3D弹性刚度
        Real lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        Real mu = E / (2.0 * (1.0 + nu));
        
        DenseMatrix D(6, 6);
        D.zero();
        
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                D(i, j) = lambda;
            }
            D(i, i) += 2.0 * mu;
        }
        
        D(3, 3) = mu;
        D(4, 4) = mu;
        D(5, 5) = mu;
        
        return D;
    } else {
        // 2D平面应变
        Real factor = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
        DenseMatrix D(3, 3);
        D.zero();
        
        D(0, 0) = factor * (1.0 - nu);
        D(0, 1) = factor * nu;
        D(1, 0) = factor * nu;
        D(1, 1) = factor * (1.0 - nu);
        D(2, 2) = factor * (1.0 - 2.0 * nu) / 2.0;
        
        return D;
    }
}

DenseMatrix J2Plasticity::consistentTangent(
    const Vector& stress,
    Real delta_gamma,
    Real equiv_plastic_strain
) const {
    // 一致性切线刚度（简化版本）
    // 完整实现需要复杂的张量运算
    // 这里返回割线刚度作为近似
    
    Real E = getParameter("E");
    Real nu = getParameter("nu");
    Real G = E / (2.0 * (1.0 + nu));
    Real H = getParameter("H");
    
    // 割线模量
    Real E_t = E * H / (E + H);  // 切线模量近似
    
    // 返回修正的弹性刚度
    DenseMatrix D_ep = elasticTensor();
    
    // 修正因子（简化）
    Real factor = 1.0 / (1.0 + 3.0 * G * delta_gamma / yieldStress(equiv_plastic_strain));
    
    for (std::size_t i = 0; i < strain_size_; ++i) {
        for (std::size_t j = 0; j < strain_size_; ++j) {
            D_ep(i, j) *= factor;
        }
    }
    
    return D_ep;
}

DenseMatrix J2Plasticity::consistentTangentFull(
    const Vector& stress,
    Real delta_gamma,
    Real equiv_plastic_strain
) const {
    // 完整一致性切线刚度（精确计算）
    // 参考：Simo & Hughes (1998), de Souza Neto et al. (2008)
    
    Real E = getParameter("E");
    Real nu = getParameter("nu");
    Real G = E / (2.0 * (1.0 + nu));
    Real K = E / (3.0 * (1.0 - 2.0 * nu));  // 体积模量
    Real H = getParameter("H");
    
    // 弹性刚度矩阵
    DenseMatrix D_e = elasticTensor();
    
    if (delta_gamma < 1e-12) {
        // 弹性状态
        return D_e;
    }
    
    // 计算偏应力
    Vector s = deviatoricStress(stress);
    Real p = hydrostaticPressure(stress);
    
    // von Mises 应力
    Real q = vonMisesStress(stress);
    
    if (q < 1e-10) {
        // 静水压力状态
        return D_e;
    }
    
    // 流动方向：n = √(3/2) * s / ||s||
    Vector n(strain_size_);
    Real sqrt_3_2 = std::sqrt(1.5);
    for (std::size_t i = 0; i < strain_size_; ++i) {
        n[i] = sqrt_3_2 * s[i] / q;
    }
    
    // 硬化模量（切线）
    Real h = 3.0 * G + H;
    
    // 一致性切线：D^ep = D^e - (3G)^2 / h * (n ⊗ n)
    DenseMatrix D_ep = D_e;
    
    Real factor = 9.0 * G * G / h;  // (3G)^2 / h
    
    for (std::size_t i = 0; i < strain_size_; ++i) {
        for (std::size_t j = 0; j < strain_size_; ++j) {
            D_ep(i, j) -= factor * n[i] * n[j];
        }
    }
    
    return D_ep;
}

}  // namespace constitutive
}  // namespace fem
