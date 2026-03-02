/**
 * kinematics.h - 几何非线性运动学工具
 * 
 * 提供大变形分析所需的运动学量计算：
 * - 变形梯度 F
 * - Green-Lagrange 应变 E
 * - Right Cauchy-Green 张量 C
 * - 应力转换
 */

#pragma once

#include "core/types.h"
#include "math/dense_matrix.h"
#include "math/vector.h"
#include <array>

namespace fem {

/**
 * Kinematics - 几何非线性运动学工具类
 * 
 * 提供 Total Lagrangian 公式所需的运动学量：
 * 
 * 1. 变形梯度 F = ∂x/∂X = I + ∇u
 * 2. Right Cauchy-Green 张量 C = F^T F
 * 3. Green-Lagrange 应变 E = 1/2(C - I)
 * 4. Jacobian J = det(F)
 * 5. 应力转换：S ↔ σ
 * 
 * 使用方法：
 * ```cpp
 * // 从位移梯度计算变形梯度
 * DenseMatrix grad_u = ...;  // ∇u (3x3)
 * DenseMatrix F = Kinematics::deformationGradient(grad_u);
 * 
 * // 计算 Green-Lagrange 应变
 * DenseMatrix E = Kinematics::greenLagrangeStrain(F);
 * 
 * // 应力转换：S → σ
 * DenseMatrix S = ...;  // 2nd Piola-Kirchhoff
 * DenseMatrix sigma = Kinematics::pushForwardStress(F, S);
 * ```
 */
class Kinematics {
public:
    // ═══ 变形梯度 ═══
    
    /**
     * 计算变形梯度 F = I + ∇u
     * 
     * @param grad_u 位移梯度 ∇u（3x3 矩阵）
     * @return F 变形梯度（3x3 矩阵）
     */
    static DenseMatrix deformationGradient(const DenseMatrix& grad_u);
    
    /**
     * 从位移梯度 Voigt 向量计算变形梯度
     * 
     * 3D Voigt: [∂u/∂x, ∂v/∂y, ∂w/∂z, ∂u/∂y, ∂v/∂x, ∂v/∂z, ∂w/∂y, ∂w/∂x, ∂u/∂z]
     * 2D Voigt: [∂u/∂x, ∂v/∂y, ∂u/∂y, ∂v/∂x]
     * 
     * @param grad_u_voigt 位移梯度（Voigt 向量）
     * @param dimension 维度（2D/3D）
     * @return F 变形梯度（3x3）
     */
    static DenseMatrix deformationGradientFromVoigt(const Vector& grad_u_voigt, int dimension);
    
    // ═══ 应变度量 ═══
    
    /**
     * 计算 Right Cauchy-Green 张量 C = F^T F
     * 
     * @param F 变形梯度（3x3）
     * @return C Right Cauchy-Green 张量（3x3）
     */
    static DenseMatrix rightCauchyGreen(const DenseMatrix& F);
    
    /**
     * 计算 Green-Lagrange 应变 E = 1/2(C - I) = 1/2(F^T F - I)
     * 
     * @param F 变形梯度（3x3）
     * @return E Green-Lagrange 应变（3x3）
     */
    static DenseMatrix greenLagrangeStrain(const DenseMatrix& F);
    
    /**
     * Green-Lagrange 应变转 Voigt 向量
     * 
     * 3D: [E_xx, E_yy, E_zz, 2E_xy, 2E_yz, 2E_xz]
     * 2D: [E_xx, E_yy, E_zz, 2E_xy]
     * 
     * @param E Green-Lagrange 应变张量（3x3）
     * @param dimension 维度
     * @return E_voigt 应变向量（Voigt 记号，工程应变）
     */
    static Vector strainToVoigt(const DenseMatrix& E, int dimension);
    
    /**
     * Voigt 向量转 Green-Lagrange 应变张量
     * 
     * @param E_voigt 应变向量（Voigt 记号）
     * @param dimension 维度
     * @return E 应变张量（3x3）
     */
    static DenseMatrix voigtToStrain(const Vector& E_voigt, int dimension);
    
    // ═══ Jacobian ═══
    
    /**
     * 计算 Jacobian J = det(F)
     * 
     * @param F 变形梯度（3x3）
     * @return J 体积比
     */
    static Real jacobian(const DenseMatrix& F);
    
    // ═══ 应力转换 ═══
    
    /**
     * Push-forward：2nd Piola-Kirchhoff → Cauchy 应力
     * 
     * σ = J^(-1) F S F^T
     * 
     * @param F 变形梯度（3x3）
     * @param S 2nd Piola-Kirchhoff 应力（3x3）
     * @return σ Cauchy 应力（3x3）
     */
    static DenseMatrix pushForwardStress(const DenseMatrix& F, const DenseMatrix& S);
    
    /**
     * Pull-back：Cauchy → 2nd Piola-Kirchhoff 应力
     * 
     * S = J F^(-1) σ F^(-T)
     * 
     * @param F 变形梯度（3x3）
     * @param sigma Cauchy 应力（3x3）
     * @return S 2nd Piola-Kirchhoff 应力（3x3）
     */
    static DenseMatrix pullBackStress(const DenseMatrix& F, const DenseMatrix& sigma);
    
    /**
     * 应力张量转 Voigt 向量
     * 
     * 3D: [σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_xz]
     * 2D: [σ_xx, σ_yy, σ_zz, σ_xy]
     * 
     * @param stress 应力张量（3x3）
     * @param dimension 维度
     * @return stress_voigt 应力向量
     */
    static Vector stressToVoigt(const DenseMatrix& stress, int dimension);
    
    /**
     * Voigt 向量转应力张量
     * 
     * @param stress_voigt 应力向量
     * @param dimension 维度
     * @return stress 应力张量（3x3）
     */
    static DenseMatrix voigtToStress(const Vector& stress_voigt, int dimension);
    
    // ═══ 不变量 ═══
    
    /**
     * 第一不变量 I_1 = tr(C)
     */
    static Real firstInvariant(const DenseMatrix& C);
    
    /**
     * 第二不变量 I_2 = 1/2[(tr C)^2 - tr(C^2)]
     */
    static Real secondInvariant(const DenseMatrix& C);
    
    /**
     * 第三不变量 I_3 = det(C)
     */
    static Real thirdInvariant(const DenseMatrix& C);
    
    // ═══ 辅助函数 ═══
    
    /**
     * 3x3 矩阵求逆
     */
    static DenseMatrix inverse3x3(const DenseMatrix& A);
    
    /**
     * 3x3 矩阵行列式
     */
    static Real det3x3(const DenseMatrix& A);
};

}  // namespace fem
