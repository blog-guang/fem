#pragma once

#include "material/material.h"

namespace fem {
namespace constitutive {

/**
 * IsotropicElastic: 各向同性线弹性材料
 * 
 * 胡克定律：σ = D : ε
 * 
 * 参数：
 * - E  : 杨氏模量 (Young's modulus)
 * - nu : 泊松比 (Poisson's ratio)
 * 
 * 约束：
 * - E > 0
 * - -1 < nu < 0.5 (3D)
 * - -1 < nu < 1 (2D平面应力)
 */
class IsotropicElastic : public Material {
public:
    // ═══ 构造函数 ═══
    
    /**
     * 构造各向同性弹性材料
     * 
     * @param E         杨氏模量
     * @param nu        泊松比
     * @param dimension 维度：2 (平面应力/应变), 3 (三维)
     * @param plane_stress 2D情况下是否为平面应力（默认true）
     */
    IsotropicElastic(Real E, Real nu, 
                     int dimension = 3, 
                     bool plane_stress = true);
    
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
    
    // ═══ 参数验证 ═══
    void validateParameters() const override;
    
    // ═══ 工具函数 ═══
    
    /**
     * 获取拉梅常数
     */
    Real lambda() const;  // 第一拉梅常数
    Real mu() const;      // 剪切模量（第二拉梅常数）
    
    /**
     * 获取体积模量
     */
    Real bulkModulus() const;
    
    /**
     * 获取弹性刚度矩阵（常数，与应变无关）
     */
    DenseMatrix elasticityTensor() const;
    
private:
    int dimension_;       // 2D or 3D
    bool plane_stress_;   // 平面应力/应变标志
    
    // 计算弹性刚度矩阵（内部辅助函数）
    DenseMatrix buildElasticityTensor() const;
};

}  // namespace constitutive
}  // namespace fem
