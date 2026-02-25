#pragma once

#include "material/material.h"

namespace fem {
namespace constitutive {

/**
 * J2PlasticityKinematic: J2塑性理论（运动硬化）
 * 
 * 特性：
 * - von Mises屈服准则
 * - 等向硬化 + 运动硬化（混合硬化）
 * - Armstrong-Frederick运动硬化规则
 * - 包辛格效应（Bauschinger effect）
 * 
 * 参数：
 * - E         : 杨氏模量
 * - nu        : 泊松比
 * - sigma_y0  : 初始屈服应力
 * - H_iso     : 等向硬化模量
 * - H_kin     : 运动硬化模量
 * - beta      : 运动硬化恢复参数（0=线性，>0=非线性）
 * 
 * 屈服函数：
 *   f = √(3/2 (s-α):(s-α)) - (σ_y0 + H_iso * ε_p^eq)
 *   其中 α 是背应力（运动硬化变量）
 * 
 * 硬化规则：
 *   dα/dε_p = (2/3) H_kin * dε_p - beta * α * d(ε_p^eq)
 *   (Armstrong-Frederick模型)
 */
class J2PlasticityKinematic : public Material {
public:
    // ═══ 构造函数 ═══
    
    /**
     * 构造运动硬化J2塑性材料
     * 
     * @param E         杨氏模量
     * @param nu        泊松比
     * @param sigma_y0  初始屈服应力
     * @param H_iso     等向硬化模量（0=无等向硬化）
     * @param H_kin     运动硬化模量（0=无运动硬化）
     * @param beta      恢复参数（0=线性，>0=非线性）
     * @param dimension 2D/3D
     */
    J2PlasticityKinematic(Real E, Real nu, 
                         Real sigma_y0, 
                         Real H_iso = 0.0,
                         Real H_kin = 0.0,
                         Real beta = 0.0,
                         int dimension = 3);
    
    // ═══ 核心接口实现 ═══
    
    void computeStress(
        const Vector& strain_inc, 
        Vector& stress, 
        StateVariables& state
    ) override;
    
    void computeTangent(
        const Vector& strain,
        DenseMatrix& D_mat,
        const StateVariables& state
    ) override;
    
    Real strainEnergy(
        const Vector& strain,
        const StateVariables& state
    ) const override;
    
    StateVariables createState() const override;
    
    std::string typeName() const override;
    
    void validateParameters() const override;
    
    // ═══ 塑性算法辅助函数 ═══
    
    /**
     * von Mises等效应力（考虑背应力）
     * q = √(3/2 * (s-α):(s-α))
     */
    Real vonMisesStress(const Vector& stress, const Vector& back_stress) const;
    
    /**
     * 屈服函数值
     * f = q - σ_y(ε_p^eq)
     */
    Real yieldFunction(Real equiv_stress, Real equiv_plastic_strain) const;
    
    /**
     * 当前屈服应力（仅等向硬化部分）
     * σ_y = σ_y0 + H_iso * ε_p^eq
     */
    Real yieldStress(Real equiv_plastic_strain) const;

private:
    int dimension_;
    
    // ═══ 内部辅助函数 ═══
    
    /**
     * 返回映射算法（含运动硬化）
     */
    void returnMapping(
        Vector& stress_trial,
        StateVariables& state,
        DenseMatrix& D_ep
    );
    
    /**
     * 计算偏应力张量
     */
    Vector deviatoricStress(const Vector& stress) const;
    
    /**
     * 计算静水压力
     */
    Real hydrostaticPressure(const Vector& stress) const;
    
    /**
     * 计算弹性刚度矩阵
     */
    DenseMatrix elasticTensor() const;
    
    /**
     * 计算一致性弹塑性切线刚度（运动硬化）
     */
    DenseMatrix consistentTangent(
        const Vector& stress,
        const Vector& back_stress,
        Real delta_gamma,
        Real equiv_plastic_strain
    ) const;
};

}  // namespace constitutive
}  // namespace fem
