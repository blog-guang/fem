/**
 * neo_hookean.cpp - Neo-Hookean 超弹性材料实现
 */

#include "material/neo_hookean.h"
#include "core/kinematics.h"
#include "core/logger.h"
#include <cmath>
#include <stdexcept>

namespace fem {
namespace constitutive {

// ═══════════════════════════════════════════════════════════
// 构造函数
// ═══════════════════════════════════════════════════════════

NeoHookean::NeoHookean(Real C10, Real D1, int dimension)
    : C10_(C10), D1_(D1), dimension_(dimension)
{
    validateParameters();
}

NeoHookean::NeoHookean(Real E, Real nu, int dimension, bool use_engineering_params)
    : dimension_(dimension)
{
    if (!use_engineering_params) {
        throw std::invalid_argument("use_engineering_params must be true");
    }
    
    computeParametersFromEngineering(E, nu);
    validateParameters();
}

void NeoHookean::computeParametersFromEngineering(Real E, Real nu) {
    // C₁₀ = E / (4(1 + ν))
    C10_ = E / (4.0 * (1.0 + nu));
    
    // D₁ = 6(1 - 2ν) / E
    // 注意：nu → 0.5 时 D₁ → 0（不可压）
    if (nu >= 0.5) {
        throw std::invalid_argument(
            "Neo-Hookean: nu must be < 0.5 (compressible material)");
    }
    D1_ = 6.0 * (1.0 - 2.0 * nu) / E;
}

// ═══════════════════════════════════════════════════════════
// 核心接口实现
// ═══════════════════════════════════════════════════════════

void NeoHookean::computeStress(
    const Vector& strain_inc,
    Vector& stress,
    StateVariables& state)
{
    // Neo-Hookean 应力计算
    // 
    // 输入 strain_inc:
    //   - 小变形模式：线性应变 ε
    //   - 大变形模式：Green-Lagrange 应变 E
    // 
    // 策略：自动识别应变大小
    //   - |ε| < 0.01: 使用线弹性近似
    //   - |ε| >= 0.01: 使用大变形公式
    
    // 计算应变范数
    Real strain_norm = 0.0;
    for (std::size_t i = 0; i < strain_inc.size(); i++) {
        strain_norm += strain_inc[i] * strain_inc[i];
    }
    strain_norm = std::sqrt(strain_norm);
    
    // 计算等效的 Lamé 参数
    Real G = 2.0 * C10_;  // 剪切模量
    Real K = 2.0 / D1_;   // 体积模量
    Real lambda = K - (2.0/3.0) * G;  // Lamé 第一参数
    
    if (strain_norm < 0.01) {
        // 小变形：线弹性近似 σ = D : ε
        compute_stress_small_strain(strain_inc, G, lambda, stress);
    } else {
        // 大变形：Neo-Hookean 超弹性
        // 从 Green-Lagrange 应变计算 Cauchy 应力
        compute_stress_large_strain(strain_inc, stress, state);
    }
    
    // 更新状态变量（用于能量计算）
    state.setScalar("W", strainEnergy(strain_inc, state));
}

void NeoHookean::compute_stress_small_strain(
    const Vector& strain,
    Real G,
    Real lambda,
    Vector& stress) const
{
    // 线弹性本构关系
    if (dimension_ == 3) {
        // 3D 应力应变关系
        Real trace_eps = strain[0] + strain[1] + strain[2];
        
        stress.resize(6);
        stress[0] = lambda * trace_eps + 2.0 * G * strain[0];  // σ_xx
        stress[1] = lambda * trace_eps + 2.0 * G * strain[1];  // σ_yy
        stress[2] = lambda * trace_eps + 2.0 * G * strain[2];  // σ_zz
        stress[3] = G * strain[3];  // σ_xy
        stress[4] = G * strain[4];  // σ_yz
        stress[5] = G * strain[5];  // σ_xz
    } else {
        // 2D 平面应变
        Real trace_eps = strain[0] + strain[1];
        
        stress.resize(4);
        stress[0] = lambda * trace_eps + 2.0 * G * strain[0];  // σ_xx
        stress[1] = lambda * trace_eps + 2.0 * G * strain[1];  // σ_yy
        stress[2] = lambda * trace_eps;  // σ_zz（平面应变约束）
        stress[3] = G * strain[2];  // σ_xy
    }
}

void NeoHookean::compute_stress_large_strain(
    const Vector& E_voigt,
    Vector& stress,
    StateVariables& state) const
{
    // 从 Green-Lagrange 应变计算 Cauchy 应力
    // 
    // 步骤：
    // 1. E (Voigt) → E (tensor)
    // 2. E → C = 2E + I
    // 3. C → F (通过极分解或直接计算)
    // 4. F → I₁, J
    // 5. 计算 S (2nd Piola-Kirchhoff)
    // 6. S → σ (Cauchy)
    
    // 简化实现：对于小到中等应变，使用线弹性近似
    // 完整实现需要存储变形梯度 F
    
    Real G = 2.0 * C10_;
    Real K = 2.0 / D1_;
    Real lambda = K - (2.0/3.0) * G;
    
    compute_stress_small_strain(E_voigt, G, lambda, stress);
    
    // TODO: 完整的大变形实现
    // DenseMatrix E = Kinematics::voigtToStrain(E_voigt, dimension_);
    // DenseMatrix C = 2 * E + I;
    // Real I1 = firstInvariant(C);
    // Real J = sqrt(thirdInvariant(C));
    // DenseMatrix S = compute2ndPiolaKirchhoff_from_C(C);
    // DenseMatrix F = compute_F_from_C(C);  // 需要极分解
    // DenseMatrix sigma = pushForwardStress(F, S);
    // stress = stressToVoigt(sigma, dimension_);
}

void NeoHookean::computeTangent(
    const Vector& strain,
    DenseMatrix& D_mat,
    const StateVariables& state)
{
    // 小变形近似：使用线弹性切线刚度
    
    Real G = 2.0 * C10_;
    Real K = 2.0 / D1_;
    Real lambda = K - (2.0/3.0) * G;
    
    if (dimension_ == 3) {
        // 3D 刚度矩阵（6x6）
        D_mat.resize(6, 6);
        D_mat.fill(0.0);
        
        Real C11 = lambda + 2.0 * G;
        Real C12 = lambda;
        
        // 对角块（法向分量）
        D_mat(0, 0) = C11; D_mat(0, 1) = C12; D_mat(0, 2) = C12;
        D_mat(1, 0) = C12; D_mat(1, 1) = C11; D_mat(1, 2) = C12;
        D_mat(2, 0) = C12; D_mat(2, 1) = C12; D_mat(2, 2) = C11;
        
        // 剪切分量
        D_mat(3, 3) = G;
        D_mat(4, 4) = G;
        D_mat(5, 5) = G;
    } else {
        // 2D 平面应变（4x4）
        D_mat.resize(4, 4);
        D_mat.fill(0.0);
        
        Real C11 = lambda + 2.0 * G;
        Real C12 = lambda;
        
        D_mat(0, 0) = C11; D_mat(0, 1) = C12;
        D_mat(1, 0) = C12; D_mat(1, 1) = C11;
        D_mat(2, 2) = C11;  // σ_zz 分量（平面应变）
        D_mat(3, 3) = G;
    }
}

Real NeoHookean::strainEnergy(
    const Vector& strain,
    const StateVariables& state) const
{
    // 小变形近似：弹性应变能
    // W = (1/2) ε : D : ε
    
    Real G = 2.0 * C10_;
    Real K = 2.0 / D1_;
    
    Real W = 0.0;
    
    if (dimension_ == 3) {
        Real eps_vol = strain[0] + strain[1] + strain[2];  // 体积应变
        Real eps_dev[3] = {
            strain[0] - eps_vol/3.0,
            strain[1] - eps_vol/3.0,
            strain[2] - eps_vol/3.0
        };
        
        // W = K/2 * eps_vol^2 + G * (eps_dev : eps_dev + gamma : gamma)
        W = 0.5 * K * eps_vol * eps_vol;
        W += G * (eps_dev[0]*eps_dev[0] + eps_dev[1]*eps_dev[1] + eps_dev[2]*eps_dev[2]);
        W += 0.5 * G * (strain[3]*strain[3] + strain[4]*strain[4] + strain[5]*strain[5]);
    } else {
        Real eps_vol = strain[0] + strain[1];
        Real eps_dev[2] = {
            strain[0] - eps_vol/2.0,
            strain[1] - eps_vol/2.0
        };
        
        W = 0.5 * K * eps_vol * eps_vol;
        W += G * (eps_dev[0]*eps_dev[0] + eps_dev[1]*eps_dev[1]);
        W += 0.5 * G * strain[2] * strain[2];
    }
    
    return W;
}

StateVariables NeoHookean::createState() const {
    StateVariables state;
    state.setScalar("W", 0.0);  // 应变能
    return state;
}

std::string NeoHookean::typeName() const {
    return "NeoHookean_" + std::to_string(dimension_) + "D";
}

void NeoHookean::validateParameters() const {
    if (C10_ <= 0.0) {
        throw std::invalid_argument("Neo-Hookean: C10 must be > 0");
    }
    if (D1_ <= 0.0) {
        throw std::invalid_argument("Neo-Hookean: D1 must be > 0");
    }
    if (dimension_ != 2 && dimension_ != 3) {
        throw std::invalid_argument("Neo-Hookean: dimension must be 2 or 3");
    }
}

// ═══════════════════════════════════════════════════════════
// 超弹性辅助函数（大变形版本 - 待完善）
// ═══════════════════════════════════════════════════════════

DenseMatrix NeoHookean::computeRightCauchyGreen(const DenseMatrix& F) {
    // C = F^T F
    DenseMatrix FT = F.transpose();
    return FT * F;
}

Real NeoHookean::firstInvariant(const DenseMatrix& C) {
    // I₁ = tr(C) = C_11 + C_22 + C_33
    return C(0, 0) + C(1, 1) + C(2, 2);
}

Real NeoHookean::jacobian(const DenseMatrix& F) {
    // J = det(F)
    // 3x3 行列式
    if (F.rows() != 3 || F.cols() != 3) {
        throw std::invalid_argument("jacobian: F must be 3x3");
    }
    
    Real det = F(0,0) * (F(1,1)*F(2,2) - F(1,2)*F(2,1))
             - F(0,1) * (F(1,0)*F(2,2) - F(1,2)*F(2,0))
             + F(0,2) * (F(1,0)*F(2,1) - F(1,1)*F(2,0));
    
    return det;
}

DenseMatrix NeoHookean::compute2ndPiolaKirchhoff(const DenseMatrix& F) const {
    // 计算 2nd Piola-Kirchhoff 应力 S = 2 ∂W/∂C
    // 
    // Neo-Hookean 应变能：
    //   W = C₁₀(Ī₁ - 3) + (1/D₁)(J - 1)²
    // 
    // 其中：
    //   Ī₁ = J^(-2/3) I₁
    //   I₁ = tr(C)
    //   J = det(F)
    //   C = F^T F
    
    // 计算 C 和不变量
    DenseMatrix C = computeRightCauchyGreen(F);
    Real I1 = firstInvariant(C);
    Real J = jacobian(F);
    
    if (J <= 0.0) {
        throw std::runtime_error("Invalid Jacobian (J <= 0): element inversion");
    }
    
    // 修正的第一不变量
    Real I1_bar = std::pow(J, -2.0/3.0) * I1;
    
    // ∂W/∂C 的计算
    // 
    // ∂W/∂I₁ = C₁₀ * J^(-2/3)
    // ∂W/∂J = 2/D₁ * (J - 1)
    // 
    // S = 2 ∂W/∂C
    //   = 2 * (∂W/∂I₁ * ∂I₁/∂C + ∂W/∂I₁_bar * ∂I₁_bar/∂C + ∂W/∂J * ∂J/∂C)
    // 
    // 简化（不可压 Neo-Hookean）：
    //   S = 2 C₁₀ J^(-2/3) (I - 1/3 I₁ C^(-1)) + 2/D₁ (J - 1) J C^(-1)
    
    // 计算 C^(-1)
    DenseMatrix C_inv = Kinematics::inverse3x3(C);
    
    // 单位矩阵
    DenseMatrix I_mat(3, 3);
    I_mat.fill(0.0);
    I_mat(0, 0) = 1.0;
    I_mat(1, 1) = 1.0;
    I_mat(2, 2) = 1.0;
    
    // 偏应力部分（等容）
    Real factor1 = 2.0 * C10_ * std::pow(J, -2.0/3.0);
    DenseMatrix S_dev = I_mat - C_inv * (I1 / 3.0);
    S_dev = S_dev * factor1;
    
    // 体积部分
    Real factor2 = 2.0 / D1_ * (J - 1.0) * J;
    DenseMatrix S_vol = C_inv * factor2;
    
    // 总的 2nd Piola-Kirchhoff 应力
    DenseMatrix S = S_dev + S_vol;
    
    return S;
}

DenseMatrix NeoHookean::computeCauchyStress(const DenseMatrix& F, const DenseMatrix& S) {
    // σ = J⁻¹ F S F^T
    Real J = jacobian(F);
    DenseMatrix FSF = F * S * F.transpose();
    return FSF * (1.0 / J);
}

void NeoHookean::computeStressFromF(const DenseMatrix& F, Vector& stress_voigt) const {
    // 完整的大变形应力计算
    // 1. F → S (2nd Piola-Kirchhoff)
    // 2. S → σ (Cauchy)
    // 3. σ → Voigt
    
    DenseMatrix S = compute2ndPiolaKirchhoff(F);
    DenseMatrix sigma = computeCauchyStress(F, S);
    stress_voigt = tensorToVoigt(sigma);
}

Vector NeoHookean::tensorToVoigt(const DenseMatrix& stress_tensor) const {
    Vector stress_voigt;
    
    if (dimension_ == 3) {
        stress_voigt.resize(6);
        stress_voigt[0] = stress_tensor(0, 0);  // σ_xx
        stress_voigt[1] = stress_tensor(1, 1);  // σ_yy
        stress_voigt[2] = stress_tensor(2, 2);  // σ_zz
        stress_voigt[3] = stress_tensor(0, 1);  // σ_xy
        stress_voigt[4] = stress_tensor(1, 2);  // σ_yz
        stress_voigt[5] = stress_tensor(0, 2);  // σ_xz
    } else {
        stress_voigt.resize(4);
        stress_voigt[0] = stress_tensor(0, 0);  // σ_xx
        stress_voigt[1] = stress_tensor(1, 1);  // σ_yy
        stress_voigt[2] = stress_tensor(2, 2);  // σ_zz
        stress_voigt[3] = stress_tensor(0, 1);  // σ_xy
    }
    
    return stress_voigt;
}

DenseMatrix NeoHookean::voigtToTensor(const Vector& strain_voigt) const {
    DenseMatrix strain_tensor(3, 3);
    strain_tensor.fill(0.0);
    
    if (dimension_ == 3) {
        strain_tensor(0, 0) = strain_voigt[0];
        strain_tensor(1, 1) = strain_voigt[1];
        strain_tensor(2, 2) = strain_voigt[2];
        strain_tensor(0, 1) = strain_tensor(1, 0) = strain_voigt[3] / 2.0;
        strain_tensor(1, 2) = strain_tensor(2, 1) = strain_voigt[4] / 2.0;
        strain_tensor(0, 2) = strain_tensor(2, 0) = strain_voigt[5] / 2.0;
    } else {
        strain_tensor(0, 0) = strain_voigt[0];
        strain_tensor(1, 1) = strain_voigt[1];
        strain_tensor(2, 2) = strain_voigt[2];
        strain_tensor(0, 1) = strain_tensor(1, 0) = strain_voigt[2] / 2.0;
    }
    
    return strain_tensor;
}

}  // namespace constitutive
}  // namespace fem
