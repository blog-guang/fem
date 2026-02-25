#include "j2_plasticity_kinematic.h"
#include <cmath>
#include <stdexcept>

namespace fem {
namespace constitutive {

J2PlasticityKinematic::J2PlasticityKinematic(Real E, Real nu, 
                                             Real sigma_y0, 
                                             Real H_iso,
                                             Real H_kin,
                                             Real beta,
                                             int dimension)
    : Material(dimension == 3 ? 6 : 3),
      dimension_(dimension)
{
    setParameter("E", E);
    setParameter("nu", nu);
    setParameter("sigma_y0", sigma_y0);
    setParameter("H_iso", H_iso);
    setParameter("H_kin", H_kin);
    setParameter("beta", beta);
    validateParameters();
}

void J2PlasticityKinematic::computeStress(
    const Vector& strain_inc, 
    Vector& stress, 
    StateVariables& state
) {
    // 初始化
    if (stress.size() != strain_size_) {
        stress.resize(strain_size_, 0.0);
    }
    
    // 确保背应力已初始化
    if (state.back_stress.size() != strain_size_) {
        state.back_stress.resize(strain_size_, 0.0);
    }
    
    // 1. 弹性预测
    DenseMatrix D_elastic = elasticTensor();
    Vector delta_stress(strain_size_, 0.0);
    
    for (std::size_t i = 0; i < strain_size_; ++i) {
        for (std::size_t j = 0; j < strain_size_; ++j) {
            delta_stress[i] += D_elastic(i, j) * strain_inc[j];
        }
    }
    
    Vector stress_trial = stress + delta_stress;
    
    // 2. 检查屈服（考虑背应力）
    Real q_trial = vonMisesStress(stress_trial, state.back_stress);
    Real sigma_y = yieldStress(state.equiv_plastic_strain);
    Real f_trial = q_trial - sigma_y;
    
    if (f_trial <= 0.0) {
        // 弹性步
        stress = stress_trial;
    } else {
        // 塑性步：返回映射
        DenseMatrix D_ep;
        returnMapping(stress_trial, state, D_ep);
        stress = stress_trial;
    }
}

void J2PlasticityKinematic::computeTangent(
    const Vector& /*strain*/,
    DenseMatrix& D_mat,
    const StateVariables& state
) {
    // 简化：返回弹性刚度
    // 完整实现需要根据当前状态判断是否塑性并返回一致性切线
    D_mat = elasticTensor();
}

Real J2PlasticityKinematic::strainEnergy(
    const Vector& strain,
    const StateVariables& state
) const {
    // 弹性应变能
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

StateVariables J2PlasticityKinematic::createState() const {
    StateVariables state(strain_size_);
    state.back_stress.resize(strain_size_, 0.0);
    return state;
}

std::string J2PlasticityKinematic::typeName() const {
    std::string name = "J2PlasticityKinematic";
    if (dimension_ == 3) {
        name += " (3D)";
    } else {
        name += " (2D)";
    }
    Real H_iso = getParameter("H_iso");
    Real H_kin = getParameter("H_kin");
    if (H_iso > 0.0 && H_kin > 0.0) {
        name += " [混合硬化]";
    } else if (H_kin > 0.0) {
        name += " [纯运动硬化]";
    } else {
        name += " [纯等向硬化]";
    }
    return name;
}

void J2PlasticityKinematic::validateParameters() const {
    Real E = getParameter("E");
    Real nu = getParameter("nu");
    Real sigma_y0 = getParameter("sigma_y0");
    Real H_iso = getParameter("H_iso");
    Real H_kin = getParameter("H_kin");
    Real beta = getParameter("beta");
    
    if (E <= 0.0) {
        throw std::invalid_argument("Young's modulus E must be positive");
    }
    if (nu <= -1.0 || nu >= 0.5) {
        throw std::invalid_argument("Poisson's ratio must be in (-1, 0.5)");
    }
    if (sigma_y0 <= 0.0) {
        throw std::invalid_argument("Initial yield stress sigma_y0 must be positive");
    }
    if (H_iso < 0.0 || H_kin < 0.0) {
        throw std::invalid_argument("Hardening moduli must be non-negative");
    }
    if (beta < 0.0) {
        throw std::invalid_argument("Recovery parameter beta must be non-negative");
    }
}

Real J2PlasticityKinematic::vonMisesStress(const Vector& stress, 
                                          const Vector& back_stress) const {
    // 相对偏应力：η = s - α
    Vector s = deviatoricStress(stress);
    
    // η = s - α
    Vector eta(strain_size_);
    for (std::size_t i = 0; i < strain_size_; ++i) {
        eta[i] = s[i] - back_stress[i];
    }
    
    // q = √(3/2 * η:η)
    Real eta_norm_sq = 0.0;
    if (dimension_ == 3) {
        eta_norm_sq = eta[0]*eta[0] + eta[1]*eta[1] + eta[2]*eta[2] 
                    + 2.0 * (eta[3]*eta[3] + eta[4]*eta[4] + eta[5]*eta[5]);
    } else {
        eta_norm_sq = eta[0]*eta[0] + eta[1]*eta[1] + 2.0 * eta[2]*eta[2];
    }
    
    return std::sqrt(1.5 * eta_norm_sq);
}

Real J2PlasticityKinematic::yieldFunction(Real equiv_stress, 
                                         Real equiv_plastic_strain) const {
    return equiv_stress - yieldStress(equiv_plastic_strain);
}

Real J2PlasticityKinematic::yieldStress(Real equiv_plastic_strain) const {
    Real sigma_y0 = getParameter("sigma_y0");
    Real H_iso = getParameter("H_iso");
    return sigma_y0 + H_iso * equiv_plastic_strain;
}

void J2PlasticityKinematic::returnMapping(
    Vector& stress_trial,
    StateVariables& state,
    DenseMatrix& D_ep
) {
    Real E = getParameter("E");
    Real nu = getParameter("nu");
    Real G = E / (2.0 * (1.0 + nu));
    Real H_iso = getParameter("H_iso");
    Real H_kin = getParameter("H_kin");
    Real beta = getParameter("beta");
    
    // 1. 分解试探应力
    Real p = hydrostaticPressure(stress_trial);
    Vector s_trial = deviatoricStress(stress_trial);
    
    // 相对偏应力
    Vector eta_trial(strain_size_);
    for (std::size_t i = 0; i < strain_size_; ++i) {
        eta_trial[i] = s_trial[i] - state.back_stress[i];
    }
    
    Real q_trial = vonMisesStress(stress_trial, state.back_stress);
    
    // 2. 计算塑性乘子（简化：不考虑非线性恢复）
    Real sigma_y = yieldStress(state.equiv_plastic_strain);
    Real denom = 3.0 * G + H_iso + H_kin;
    Real delta_gamma = (q_trial - sigma_y) / denom;
    
    if (delta_gamma < 0.0) {
        delta_gamma = 0.0;
    }
    
    // 3. 更新背应力（Armstrong-Frederick）
    // Δα = (2/3) H_kin * Δγ * n - beta * α * Δγ
    // 其中 n = η / ||η||
    
    Real sqrt_3_2 = std::sqrt(1.5);
    Real eta_norm = q_trial / sqrt_3_2;
    
    if (eta_norm > 1e-10) {
        for (std::size_t i = 0; i < strain_size_; ++i) {
            Real n_i = eta_trial[i] / eta_norm;
            Real delta_alpha = (2.0/3.0) * H_kin * delta_gamma * n_i
                              - beta * state.back_stress[i] * delta_gamma;
            state.back_stress[i] += delta_alpha;
        }
    }
    
    // 4. 修正应力（径向返回）
    Real factor = 1.0 - (3.0 * G * delta_gamma) / q_trial;
    if (factor < 0.0) factor = 0.0;
    
    Vector s_corrected(strain_size_);
    for (std::size_t i = 0; i < strain_size_; ++i) {
        // s_new = α + (1-3G*Δγ/q) * (s_trial - α)
        s_corrected[i] = state.back_stress[i] + factor * eta_trial[i];
    }
    
    // 5. 重构应力
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
    
    // 6. 更新塑性应变
    Real q_new = vonMisesStress(stress_trial, state.back_stress);
    Real norm_eta = q_new / sqrt_3_2;
    
    if (norm_eta > 1e-10) {
        Vector eta_new(strain_size_);
        Vector s_new = deviatoricStress(stress_trial);
        for (std::size_t i = 0; i < strain_size_; ++i) {
            eta_new[i] = s_new[i] - state.back_stress[i];
        }
        
        for (std::size_t i = 0; i < strain_size_; ++i) {
            Real flow_dir = eta_new[i] / norm_eta;
            state.plastic_strain[i] += sqrt_3_2 * delta_gamma * flow_dir;
        }
    }
    
    state.equiv_plastic_strain += delta_gamma;
    
    // 7. 一致性切线（简化）
    D_ep = consistentTangent(stress_trial, state.back_stress, 
                            delta_gamma, state.equiv_plastic_strain);
}

Vector J2PlasticityKinematic::deviatoricStress(const Vector& stress) const {
    Real p = hydrostaticPressure(stress);
    
    Vector s = stress;
    if (dimension_ == 3) {
        s[0] -= p;
        s[1] -= p;
        s[2] -= p;
    } else {
        s[0] -= p;
        s[1] -= p;
    }
    
    return s;
}

Real J2PlasticityKinematic::hydrostaticPressure(const Vector& stress) const {
    if (dimension_ == 3) {
        return (stress[0] + stress[1] + stress[2]) / 3.0;
    } else {
        return (stress[0] + stress[1]) / 2.0;
    }
}

DenseMatrix J2PlasticityKinematic::elasticTensor() const {
    Real E = getParameter("E");
    Real nu = getParameter("nu");
    
    if (dimension_ == 3) {
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

DenseMatrix J2PlasticityKinematic::consistentTangent(
    const Vector& stress,
    const Vector& back_stress,
    Real delta_gamma,
    Real equiv_plastic_strain
) const {
    // 简化实现：割线刚度
    Real E = getParameter("E");
    Real nu = getParameter("nu");
    Real G = E / (2.0 * (1.0 + nu));
    Real H_iso = getParameter("H_iso");
    Real H_kin = getParameter("H_kin");
    
    DenseMatrix D_ep = elasticTensor();
    
    Real factor = 1.0 / (1.0 + (3.0*G + H_iso + H_kin) * delta_gamma / 
                                yieldStress(equiv_plastic_strain));
    
    for (std::size_t i = 0; i < strain_size_; ++i) {
        for (std::size_t j = 0; j < strain_size_; ++j) {
            D_ep(i, j) *= factor;
        }
    }
    
    return D_ep;
}

}  // namespace constitutive
}  // namespace fem
