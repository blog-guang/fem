/**
 * neo_hookean.h - Neo-Hookean 超弹性材料模型
 * 
 * 用于大变形超弹性材料（橡胶、软组织、泡沫等）
 * 
 * 理论：
 * - 基于应变能密度函数
 * - 支持不可压/近似不可压
 * - Total Lagrangian 公式
 */

#pragma once

#include "material/material.h"
#include "math/dense_matrix.h"

namespace fem {
namespace constitutive {

/**
 * NeoHookean: Neo-Hookean 超弹性模型
 * 
 * 应变能密度函数：
 *   W = C₁₀(Ī₁ - 3) + (1/D₁)(J - 1)²
 * 
 * 其中：
 *   Ī₁ = J^(-2/3) I₁  (修正的第一不变量，用于不可压性)
 *   I₁ = tr(C) = λ₁² + λ₂² + λ₃²  (右 Cauchy-Green 第一不变量)
 *   J = det(F)  (体积比)
 *   C₁₀ = G/2  (剪切模量参数)
 *   D₁ = 2/κ  (体积模量参数)
 * 
 * 参数输入：
 *   - C10: 材料参数 C₁₀ = G/2
 *   - D1:  材料参数 D₁ = 2/κ
 *   或者
 *   - E:   杨氏模量
 *   - nu:  泊松比
 *   （通过 C₁₀ = E/(4(1+ν)), D₁ = 6(1-2ν)/E 转换）
 * 
 * 用途：
 *   - 橡胶密封件
 *   - 生物软组织
 *   - 泡沫材料
 *   - 高分子聚合物
 * 
 * 限制：
 *   - 仅适用于单调加载路径
 *   - 不考虑粘弹性
 *   - 不考虑损伤
 */
class NeoHookean : public Material {
public:
    // ═══ 构造函数 ═══
    
    /**
     * 通过材料参数 C₁₀, D₁ 构造
     * 
     * @param C10 剪切参数 C₁₀ = G/2
     * @param D1  体积参数 D₁ = 2/κ
     * @param dimension 维度（2D 平面应变, 3D）
     */
    NeoHookean(Real C10, Real D1, int dimension = 3);
    
    /**
     * 通过工程参数 E, ν 构造
     * 
     * @param E 杨氏模量
     * @param nu 泊松比
     * @param dimension 维度（2D 平面应变, 3D）
     * @param use_engineering_params 标志（必须为 true）
     */
    NeoHookean(Real E, Real nu, int dimension, bool use_engineering_params);
    
    // ═══ 核心接口 ═══
    
    /**
     * 计算应力（基于变形梯度）
     * 
     * @param strain_inc 变形梯度增量（以 Voigt 形式存储）
     * @param stress 输出：Cauchy 应力
     * @param state 状态变量（存储 F, C 等）
     */
    void computeStress(
        const Vector& strain_inc,
        Vector& stress,
        StateVariables& state
    ) override;
    
    /**
     * 计算切线刚度矩阵
     * 
     * @param strain 当前应变（变形梯度）
     * @param D_mat 输出：材料切线刚度矩阵
     * @param state 状态变量
     */
    void computeTangent(
        const Vector& strain,
        DenseMatrix& D_mat,
        const StateVariables& state
    ) override;
    
    /**
     * 应变能密度
     */
    Real strainEnergy(
        const Vector& strain,
        const StateVariables& state
    ) const override;
    
    StateVariables createState() const override;
    
    std::string typeName() const override;
    
    void validateParameters() const override;
    
    // ═══ 超弹性辅助函数 ═══
    
    /**
     * 从变形梯度 F 计算右 Cauchy-Green 张量 C = F^T F
     * 
     * @param F 变形梯度（3x3 矩阵）
     * @return C 右 Cauchy-Green 张量（3x3 矩阵）
     */
    static DenseMatrix computeRightCauchyGreen(const DenseMatrix& F);
    
    /**
     * 计算第一不变量 I₁ = tr(C)
     */
    static Real firstInvariant(const DenseMatrix& C);
    
    /**
     * 计算 Jacobian J = det(F)
     */
    static Real jacobian(const DenseMatrix& F);
    
    /**
     * 计算 2nd Piola-Kirchhoff 应力 S = 2 ∂W/∂C
     * 
     * @param F 变形梯度（3x3）
     * @return S 2nd Piola-Kirchhoff 应力（3x3）
     */
    DenseMatrix compute2ndPiolaKirchhoff(const DenseMatrix& F) const;
    
    /**
     * 计算 Cauchy 应力 σ = J⁻¹ F S F^T
     * 
     * @param F 变形梯度（3x3）
     * @param S 2nd Piola-Kirchhoff 应力（3x3）
     * @return σ Cauchy 应力（3x3）
     */
    static DenseMatrix computeCauchyStress(const DenseMatrix& F, const DenseMatrix& S);
    
    /**
     * Voigt 记号转换：应力张量 → 应力向量
     * 
     * 3D: [σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_xz]
     * 2D 平面应变: [σ_xx, σ_yy, σ_zz, σ_xy]
     */
    Vector tensorToVoigt(const DenseMatrix& stress_tensor) const;
    
    /**
     * Voigt 记号转换：应变向量 → 应变张量
     */
    DenseMatrix voigtToTensor(const Vector& strain_voigt) const;
    
private:
    // ═══ 材料参数 ═══
    Real C10_;       ///< 剪切参数 C₁₀ = G/2
    Real D1_;        ///< 体积参数 D₁ = 2/κ
    int dimension_;  ///< 维度（2D/3D）
    
    // ═══ 内部辅助函数 ═══
    
    /**
     * 从工程参数计算材料参数
     */
    void computeParametersFromEngineering(Real E, Real nu);
};

}  // namespace constitutive
}  // namespace fem
