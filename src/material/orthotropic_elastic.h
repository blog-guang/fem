#pragma once

#include "material/material.h"

namespace fem {
namespace constitutive {

/**
 * OrthotropicElastic: 正交各向异性弹性材料
 * 
 * 典型应用：复合材料、木材、层压板
 * 
 * 材料主轴参数：
 * - E1, E2, E3  : 三个主方向的杨氏模量
 * - nu12, nu23, nu13 : 主泊松比
 * - G12, G23, G13    : 主剪切模量
 * 
 * 2D情况（平面应力）：
 * - E1, E2  : 面内杨氏模量
 * - nu12    : 面内泊松比
 * - G12     : 面内剪切模量
 * 
 * 坐标变换：
 * - theta : 主轴与全局坐标系的夹角（逆时针为正）
 */
class OrthotropicElastic : public Material {
public:
    // ═══ 构造函数 ═══
    
    /**
     * 2D正交各向异性材料
     * 
     * @param E1    主方向1的杨氏模量
     * @param E2    主方向2的杨氏模量
     * @param nu12  主泊松比
     * @param G12   主剪切模量
     * @param theta 主轴旋转角（度数）
     */
    OrthotropicElastic(Real E1, Real E2, Real nu12, Real G12, 
                      Real theta = 0.0);
    
    /**
     * 3D正交各向异性材料
     * 
     * @param E1, E2, E3  主方向杨氏模量
     * @param nu12, nu23, nu13  主泊松比
     * @param G12, G23, G13     主剪切模量
     */
    OrthotropicElastic(Real E1, Real E2, Real E3,
                      Real nu12, Real nu23, Real nu13,
                      Real G12, Real G23, Real G13);
    
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
    
    // ═══ 工具函数 ═══
    
    /**
     * 获取主轴刚度矩阵（局部坐标系）
     */
    DenseMatrix localStiffness() const;
    
    /**
     * 获取全局刚度矩阵（考虑旋转）
     */
    DenseMatrix globalStiffness() const;
    
    /**
     * 设置旋转角度
     */
    void setRotation(Real theta_deg);
    Real rotation() const { return theta_; }

private:
    int dimension_;  // 2 or 3
    Real theta_;     // 旋转角（弧度）
    
    // ═══ 内部辅助函数 ═══
    
    /**
     * 计算2D应力变换矩阵
     * T: 从主轴到全局坐标
     */
    DenseMatrix transformationMatrix2D() const;
    
    /**
     * 构建局部刚度矩阵（主轴坐标系）
     */
    DenseMatrix buildLocalStiffness2D() const;
    DenseMatrix buildLocalStiffness3D() const;
};

}  // namespace constitutive
}  // namespace fem
