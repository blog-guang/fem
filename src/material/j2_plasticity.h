#pragma once

#include "material/material.h"

namespace fem {
namespace constitutive {

/**
 * J2Plasticity: J2塑性理论（von Mises塑性）
 * 
 * 特性：
 * - von Mises屈服准则
 * - 等向硬化（线性/非线性）
 * - 关联流动法则
 * - 径向返回映射算法（Radial Return Mapping）
 * - 一致性切线刚度矩阵
 * 
 * 参数：
 * - E         : 杨氏模量
 * - nu        : 泊松比
 * - sigma_y0  : 初始屈服应力
 * - H         : 硬化模量（等向硬化）
 * - dimension : 2D/3D
 * 
 * 屈服函数：
 *   f = √(3/2 J2) - σ_y(ε_p^eq)
 *   其中 J2 = 0.5 * s:s (偏应力第二不变量)
 *   σ_y = σ_y0 + H * ε_p^eq (线性硬化)
 */
class J2Plasticity : public Material {
public:
    // ═══ 构造函数 ═══
    
    /**
     * 构造J2塑性材料
     * 
     * @param E         杨氏模量
     * @param nu        泊松比
     * @param sigma_y0  初始屈服应力
     * @param H         硬化模量（0=理想塑性）
     * @param dimension 2D/3D
     */
    J2Plasticity(Real E, Real nu, 
                 Real sigma_y0, 
                 Real H = 0.0,
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
     * von Mises等效应力
     * q = √(3/2 * s:s) = √(3 * J2)
     */
    Real vonMisesStress(const Vector& stress) const;
    
    /**
     * 屈服函数值
     * f = q - σ_y(ε_p^eq)
     */
    Real yieldFunction(Real equiv_stress, Real equiv_plastic_strain) const;
    
    /**
     * 当前屈服应力（含硬化）
     * σ_y = σ_y0 + H * ε_p^eq
     */
    Real yieldStress(Real equiv_plastic_strain) const;
    
private:
    int dimension_;
    
    // ═══ 内部辅助函数 ═══
    
    /**
     * 返回映射算法（Return Mapping）
     * 
     * 给定试探应力 σ_trial，修正到屈服面上
     * 更新状态变量（塑性应变、等效塑性应变）
     */
    void returnMapping(
        Vector& stress_trial,
        StateVariables& state,
        DenseMatrix& D_ep
    );
    
    /**
     * 计算偏应力张量
     * s = σ - (1/3) tr(σ) I
     */
    Vector deviatoricStress(const Vector& stress) const;
    
    /**
     * 计算应力迹（静水压力）
     * p = tr(σ) / 3
     */
    Real hydrostaticPressure(const Vector& stress) const;
    
    /**
     * 计算弹性刚度矩阵
     */
    DenseMatrix elasticTensor() const;
    
    /**
     * 计算一致性弹塑性切线刚度
     * (Consistent tangent for Newton-Raphson)
     */
    DenseMatrix consistentTangent(
        const Vector& stress,
        Real delta_gamma,
        Real equiv_plastic_strain
    ) const;
};

}  // namespace constitutive
}  // namespace fem
